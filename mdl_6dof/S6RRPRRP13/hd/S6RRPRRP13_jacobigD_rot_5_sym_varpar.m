% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPRRP13
% Use Code from Maple symbolic Code Generation
%
% geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5]';
%
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:53
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6RRPRRP13_jacobigD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP13_jacobigD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP13_jacobigD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP13_jacobigD_rot_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:52:49
% EndTime: 2019-02-26 21:52:49
% DurationCPUTime: 0.06s
% Computational Cost: add. (22->16), mult. (78->39), div. (0->0), fcn. (80->8), ass. (0->22)
t155 = sin(pkin(6));
t159 = sin(qJ(1));
t175 = t155 * t159;
t161 = cos(qJ(2));
t174 = t155 * t161;
t162 = cos(qJ(1));
t173 = t155 * t162;
t158 = sin(qJ(2));
t172 = t159 * t158;
t171 = t159 * t161;
t170 = t161 * t162;
t169 = t162 * t158;
t168 = qJD(1) * t155;
t160 = cos(qJ(4));
t167 = qJD(2) * t160;
t156 = cos(pkin(6));
t166 = t156 * t170 - t172;
t165 = t156 * t171 + t169;
t164 = t156 * t169 + t171;
t163 = -t156 * t172 + t170;
t157 = sin(qJ(4));
t1 = [0, t162 * t168, 0, -t164 * qJD(1) - t165 * qJD(2) (t165 * t157 + t160 * t175) * qJD(4) - t163 * t167 + (t157 * t173 - t166 * t160) * qJD(1), 0; 0, t159 * t168, 0, t163 * qJD(1) + t166 * qJD(2) (-t166 * t157 - t160 * t173) * qJD(4) - t164 * t167 + (t157 * t175 - t165 * t160) * qJD(1), 0; 0, 0, 0, qJD(2) * t174, -t155 * t158 * t167 + (t156 * t160 - t157 * t174) * qJD(4), 0;];
JgD_rot  = t1;

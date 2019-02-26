% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPPRR5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6]';
%
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:30
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6RRPPRR5_jacobigD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR5_jacobigD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR5_jacobigD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR5_jacobigD_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:30:48
% EndTime: 2019-02-26 21:30:48
% DurationCPUTime: 0.06s
% Computational Cost: add. (23->17), mult. (78->40), div. (0->0), fcn. (80->8), ass. (0->21)
t153 = sin(pkin(6));
t158 = cos(qJ(5));
t172 = t153 * t158;
t160 = cos(qJ(1));
t171 = t153 * t160;
t156 = sin(qJ(2));
t157 = sin(qJ(1));
t170 = t156 * t157;
t169 = t156 * t160;
t159 = cos(qJ(2));
t168 = t157 * t159;
t167 = t160 * t159;
t166 = qJD(1) * t153;
t155 = sin(qJ(5));
t165 = qJD(2) * t155;
t154 = cos(pkin(6));
t164 = t154 * t167 - t170;
t163 = -t154 * t168 - t169;
t162 = t154 * t169 + t168;
t161 = t154 * t170 - t167;
t1 = [0, t160 * t166, 0, 0, -t164 * qJD(1) + t161 * qJD(2) (-t157 * t153 * t155 - t161 * t158) * qJD(5) + t163 * t165 + (-t162 * t155 + t158 * t171) * qJD(1); 0, t157 * t166, 0, 0, t163 * qJD(1) - t162 * qJD(2) (t155 * t171 + t162 * t158) * qJD(5) + t164 * t165 + (-t161 * t155 + t157 * t172) * qJD(1); 0, 0, 0, 0, -t153 * qJD(2) * t156, t153 * t159 * t165 + (-t154 * t155 + t156 * t172) * qJD(5);];
JgD_rot  = t1;

% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRPRRP6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1]';
%
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:53
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6PRPRRP6_jacobigD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP6_jacobigD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRP6_jacobigD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRPRRP6_jacobigD_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:53:16
% EndTime: 2019-02-26 19:53:16
% DurationCPUTime: 0.04s
% Computational Cost: add. (12->11), mult. (48->29), div. (0->0), fcn. (50->8), ass. (0->16)
t157 = sin(pkin(6));
t162 = cos(qJ(4));
t170 = t157 * t162;
t163 = cos(qJ(2));
t169 = t157 * t163;
t159 = cos(pkin(6));
t161 = sin(qJ(2));
t168 = t159 * t161;
t167 = t159 * t163;
t166 = qJD(2) * t162;
t156 = sin(pkin(10));
t158 = cos(pkin(10));
t165 = -t156 * t161 + t158 * t167;
t164 = t156 * t167 + t158 * t161;
t160 = sin(qJ(4));
t1 = [0, 0, 0, -t164 * qJD(2) (t156 * t170 + t164 * t160) * qJD(4) - (-t156 * t168 + t158 * t163) * t166, 0; 0, 0, 0, t165 * qJD(2) (-t158 * t170 - t165 * t160) * qJD(4) - (t156 * t163 + t158 * t168) * t166, 0; 0, 0, 0, qJD(2) * t169, -t157 * t161 * t166 + (t159 * t162 - t160 * t169) * qJD(4), 0;];
JgD_rot  = t1;

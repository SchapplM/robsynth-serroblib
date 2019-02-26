% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRPRRP4
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1,theta3]';
%
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:52
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6PRPRRP4_jacobigD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP4_jacobigD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRP4_jacobigD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP4_jacobigD_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:52:01
% EndTime: 2019-02-26 19:52:01
% DurationCPUTime: 0.04s
% Computational Cost: add. (21->11), mult. (48->30), div. (0->0), fcn. (50->8), ass. (0->17)
t163 = pkin(11) + qJ(4);
t161 = sin(t163);
t165 = sin(pkin(6));
t176 = t165 * t161;
t167 = cos(pkin(6));
t168 = sin(qJ(2));
t175 = t167 * t168;
t169 = cos(qJ(2));
t174 = t167 * t169;
t173 = qJD(2) * t161;
t172 = qJD(2) * t165;
t164 = sin(pkin(10));
t166 = cos(pkin(10));
t171 = t164 * t169 + t166 * t175;
t170 = -t164 * t175 + t166 * t169;
t162 = cos(t163);
t1 = [0, 0, 0, t170 * qJD(2) (t170 * t162 + t164 * t176) * qJD(4) + (-t164 * t174 - t166 * t168) * t173, 0; 0, 0, 0, t171 * qJD(2) (t171 * t162 - t166 * t176) * qJD(4) + (-t164 * t168 + t166 * t174) * t173, 0; 0, 0, 0, t168 * t172, t169 * t161 * t172 + (t162 * t165 * t168 + t161 * t167) * qJD(4), 0;];
JgD_rot  = t1;

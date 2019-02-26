% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PPRPRR1
% Use Code from Maple symbolic Code Generation
%
% geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d5,d6,theta1,theta2,theta4]';
%
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:39
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6PPRPRR1_jacobig_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRPRR1_jacobig_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRPRR1_jacobig_rot_6_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:39:48
% EndTime: 2019-02-26 19:39:48
% DurationCPUTime: 0.07s
% Computational Cost: add. (51->24), mult. (146->54), div. (0->0), fcn. (206->14), ass. (0->33)
t153 = sin(pkin(11));
t155 = sin(pkin(6));
t170 = t153 * t155;
t160 = cos(pkin(6));
t169 = t153 * t160;
t158 = cos(pkin(11));
t168 = t155 * t158;
t159 = cos(pkin(7));
t167 = t155 * t159;
t166 = t158 * t160;
t151 = sin(pkin(13));
t156 = cos(pkin(13));
t162 = sin(qJ(3));
t164 = cos(qJ(3));
t165 = t164 * t151 + t162 * t156;
t150 = -t162 * t151 + t164 * t156;
t163 = cos(qJ(5));
t161 = sin(qJ(5));
t157 = cos(pkin(12));
t154 = sin(pkin(7));
t152 = sin(pkin(12));
t148 = -t152 * t169 + t158 * t157;
t147 = -t158 * t152 - t157 * t169;
t146 = t152 * t166 + t153 * t157;
t145 = -t153 * t152 + t157 * t166;
t144 = -t155 * t157 * t154 + t160 * t159;
t143 = t165 * t159;
t142 = t150 * t159;
t141 = t165 * t154;
t140 = t150 * t154;
t139 = -t147 * t154 + t153 * t167;
t138 = -t145 * t154 - t158 * t167;
t1 = [0, 0, t139, 0, -t140 * t170 - t147 * t142 + t148 * t165 (t141 * t170 + t147 * t143 + t148 * t150) * t161 - t139 * t163; 0, 0, t138, 0, t140 * t168 - t145 * t142 + t146 * t165 (-t141 * t168 + t145 * t143 + t146 * t150) * t161 - t138 * t163; 0, 0, t144, 0, -t160 * t140 + (-t142 * t157 + t152 * t165) * t155 (t160 * t141 + (t143 * t157 + t150 * t152) * t155) * t161 - t144 * t163;];
Jg_rot  = t1;

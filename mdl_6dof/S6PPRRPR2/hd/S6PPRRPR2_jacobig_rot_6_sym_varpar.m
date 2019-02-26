% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PPRRPR2
% Use Code from Maple symbolic Code Generation
%
% geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d6,theta1,theta2]';
%
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:41
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6PPRRPR2_jacobig_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRPR2_jacobig_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PPRRPR2_jacobig_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:40:57
% EndTime: 2019-02-26 19:40:57
% DurationCPUTime: 0.07s
% Computational Cost: add. (33->20), mult. (98->46), div. (0->0), fcn. (140->12), ass. (0->28)
t141 = sin(pkin(11));
t147 = cos(pkin(6));
t159 = t141 * t147;
t142 = sin(pkin(7));
t143 = sin(pkin(6));
t158 = t142 * t143;
t157 = t142 * t147;
t146 = cos(pkin(7));
t156 = t143 * t146;
t144 = cos(pkin(12));
t155 = t144 * t146;
t145 = cos(pkin(11));
t154 = t145 * t147;
t140 = sin(pkin(12));
t136 = -t140 * t141 + t144 * t154;
t153 = -t136 * t146 + t145 * t158;
t138 = -t140 * t145 - t144 * t159;
t152 = t138 * t146 + t141 * t158;
t151 = cos(qJ(3));
t150 = cos(qJ(4));
t149 = sin(qJ(3));
t148 = sin(qJ(4));
t139 = -t140 * t159 + t144 * t145;
t137 = t140 * t154 + t141 * t144;
t135 = -t144 * t158 + t146 * t147;
t134 = -t138 * t142 + t141 * t156;
t133 = -t136 * t142 - t145 * t156;
t1 = [0, 0, t134, t139 * t149 - t152 * t151, 0 (t139 * t151 + t152 * t149) * t150 + t134 * t148; 0, 0, t133, t137 * t149 + t153 * t151, 0 (t137 * t151 - t153 * t149) * t150 + t133 * t148; 0, 0, t135, -t151 * t157 + (t140 * t149 - t151 * t155) * t143, 0 (t149 * t157 + (t140 * t151 + t149 * t155) * t143) * t150 + t135 * t148;];
Jg_rot  = t1;

% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6PPRRRR3
% Use Code from Maple symbolic Code Generation
%
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d3,d4,d5,d6,theta1,theta2]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 09:27
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6PPRRRR3_jacobiR_rot_4_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRR3_jacobiR_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PPRRRR3_jacobiR_rot_4_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 09:27:11
% EndTime: 2019-02-22 09:27:11
% DurationCPUTime: 0.12s
% Computational Cost: add. (121->33), mult. (360->77), div. (0->0), fcn. (493->14), ass. (0->38)
t132 = sin(pkin(13));
t140 = cos(pkin(6));
t157 = t132 * t140;
t134 = sin(pkin(7));
t135 = sin(pkin(6));
t156 = t134 * t135;
t155 = t134 * t140;
t139 = cos(pkin(7));
t154 = t135 * t139;
t137 = cos(pkin(13));
t153 = t137 * t140;
t138 = cos(pkin(8));
t141 = sin(qJ(4));
t152 = t138 * t141;
t143 = cos(qJ(4));
t151 = t138 * t143;
t142 = sin(qJ(3));
t150 = t139 * t142;
t149 = t137 * t156;
t131 = sin(pkin(14));
t136 = cos(pkin(14));
t126 = -t132 * t131 + t136 * t153;
t127 = t131 * t153 + t132 * t136;
t144 = cos(qJ(3));
t117 = -t127 * t142 + (t126 * t139 - t149) * t144;
t133 = sin(pkin(8));
t148 = t117 * t138 + (-t126 * t134 - t137 * t154) * t133;
t129 = -t131 * t157 + t137 * t136;
t128 = -t137 * t131 - t136 * t157;
t145 = t128 * t139 + t132 * t156;
t119 = -t129 * t142 + t144 * t145;
t147 = t119 * t138 + (-t128 * t134 + t132 * t154) * t133;
t121 = t144 * t155 + (t136 * t139 * t144 - t131 * t142) * t135;
t146 = t121 * t138 + (-t136 * t156 + t140 * t139) * t133;
t122 = t142 * t155 + (t131 * t144 + t136 * t150) * t135;
t120 = t129 * t144 + t142 * t145;
t118 = t126 * t150 + t127 * t144 - t142 * t149;
t1 = [0, 0, t119 * t143 - t120 * t152, -t120 * t141 + t143 * t147, 0, 0; 0, 0, t117 * t143 - t118 * t152, -t118 * t141 + t143 * t148, 0, 0; 0, 0, t121 * t143 - t122 * t152, -t122 * t141 + t143 * t146, 0, 0; 0, 0, -t119 * t141 - t120 * t151, -t120 * t143 - t141 * t147, 0, 0; 0, 0, -t117 * t141 - t118 * t151, -t118 * t143 - t141 * t148, 0, 0; 0, 0, -t121 * t141 - t122 * t151, -t122 * t143 - t141 * t146, 0, 0; 0, 0, t120 * t133, 0, 0, 0; 0, 0, t118 * t133, 0, 0, 0; 0, 0, t122 * t133, 0, 0, 0;];
JR_rot  = t1;

% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRPPR4
% Use Code from Maple symbolic Code Generation
%
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta4]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 09:44
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6PRRPPR4_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR4_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPPR4_jacobiR_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 09:44:03
% EndTime: 2019-02-22 09:44:04
% DurationCPUTime: 0.14s
% Computational Cost: add. (116->40), mult. (337->91), div. (0->0), fcn. (474->12), ass. (0->45)
t144 = sin(pkin(11));
t154 = cos(qJ(3));
t165 = t144 * t154;
t146 = sin(pkin(6));
t151 = sin(qJ(3));
t164 = t146 * t151;
t163 = t146 * t154;
t155 = cos(qJ(2));
t162 = t146 * t155;
t147 = cos(pkin(11));
t161 = t147 * t154;
t149 = cos(pkin(6));
t152 = sin(qJ(2));
t160 = t149 * t152;
t159 = t149 * t155;
t158 = t154 * t155;
t150 = sin(qJ(6));
t153 = cos(qJ(6));
t157 = t144 * t153 - t147 * t150;
t156 = t144 * t150 + t147 * t153;
t148 = cos(pkin(10));
t145 = sin(pkin(10));
t142 = t149 * t151 + t152 * t163;
t141 = t149 * t154 - t152 * t164;
t140 = -t145 * t160 + t148 * t155;
t139 = t145 * t159 + t148 * t152;
t138 = t145 * t155 + t148 * t160;
t137 = t145 * t152 - t148 * t159;
t136 = (t144 * t152 + t147 * t158) * t146;
t135 = (t144 * t158 - t147 * t152) * t146;
t134 = t140 * t154 + t145 * t164;
t133 = -t140 * t151 + t145 * t163;
t132 = t138 * t154 - t148 * t164;
t131 = -t138 * t151 - t148 * t163;
t130 = t142 * t147 - t144 * t162;
t129 = t142 * t144 + t147 * t162;
t128 = -t139 * t161 + t140 * t144;
t127 = -t139 * t165 - t140 * t147;
t126 = -t137 * t161 + t138 * t144;
t125 = -t137 * t165 - t138 * t147;
t124 = t134 * t147 + t139 * t144;
t123 = t134 * t144 - t139 * t147;
t122 = t132 * t147 + t137 * t144;
t121 = t132 * t144 - t137 * t147;
t1 = [0, t127 * t150 + t128 * t153, t156 * t133, 0, 0, t123 * t153 - t124 * t150; 0, t125 * t150 + t126 * t153, t156 * t131, 0, 0, t121 * t153 - t122 * t150; 0, t135 * t150 + t136 * t153, t156 * t141, 0, 0, t129 * t153 - t130 * t150; 0, t127 * t153 - t128 * t150, t157 * t133, 0, 0, -t123 * t150 - t124 * t153; 0, t125 * t153 - t126 * t150, t157 * t131, 0, 0, -t121 * t150 - t122 * t153; 0, t135 * t153 - t136 * t150, t157 * t141, 0, 0, -t129 * t150 - t130 * t153; 0, t139 * t151, -t134, 0, 0, 0; 0, t137 * t151, -t132, 0, 0, 0; 0, -t151 * t162, -t142, 0, 0, 0;];
JR_rot  = t1;

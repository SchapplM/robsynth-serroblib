% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 3 (0=Basis) von
% S6RRPRRR14
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d2,d4,d5,d6,theta3]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-12-10 18:38
% Revision: bb42a8b95257d9bc83910d26e849f5825122f662 (2018-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function JR_rot = S6RRPRRR14_jacobiR_rot_3_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR14_jacobiR_rot_3_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRPRRR14_jacobiR_rot_3_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_3_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2018-12-10 18:38:21
% EndTime: 2018-12-10 18:38:21
% DurationCPUTime: 0.08s
% Computational Cost: add. (156->35), mult. (155->53), div. (0->0), fcn. (167->17), ass. (0->36)
t138 = sin(pkin(6));
t142 = sin(qJ(1));
t151 = t138 * t142;
t144 = cos(qJ(1));
t150 = t138 * t144;
t149 = pkin(6) - qJ(2);
t135 = pkin(6) + qJ(2);
t129 = cos(t135) / 0.2e1;
t132 = cos(t149);
t148 = t132 / 0.2e1 + t129;
t147 = sin(t135) / 0.2e1;
t146 = sin(t149);
t123 = t147 - t146 / 0.2e1;
t143 = cos(qJ(2));
t145 = -t144 * t123 - t142 * t143;
t116 = t142 * t123 - t144 * t143;
t141 = sin(qJ(2));
t140 = cos(pkin(7));
t139 = cos(pkin(14));
t137 = sin(pkin(7));
t136 = sin(pkin(14));
t134 = pkin(7) - pkin(14);
t133 = pkin(7) + pkin(14);
t131 = cos(t133);
t130 = sin(t134);
t126 = cos(t134) / 0.2e1;
t125 = sin(t133) / 0.2e1;
t124 = t129 - t132 / 0.2e1;
t122 = t147 + t146 / 0.2e1;
t121 = t126 - t131 / 0.2e1;
t120 = t126 + t131 / 0.2e1;
t119 = t125 - t130 / 0.2e1;
t118 = t125 + t130 / 0.2e1;
t114 = -t144 * t141 - t142 * t148;
t111 = t142 * t141 - t144 * t148;
t1 = [t111 * t119 + t121 * t150 + t139 * t145, t114 * t139 + t116 * t119, 0, 0, 0, 0; t114 * t119 - t116 * t139 + t121 * t151, -t111 * t139 + t119 * t145, 0, 0, 0, 0; 0, t124 * t119 + t122 * t139, 0, 0, 0, 0; t111 * t120 + t118 * t150 - t136 * t145, -t114 * t136 + t116 * t120, 0, 0, 0, 0; t114 * t120 + t116 * t136 + t118 * t151, t111 * t136 + t120 * t145, 0, 0, 0, 0; 0, t124 * t120 - t122 * t136, 0, 0, 0, 0; -t111 * t137 + t140 * t150, -t116 * t137, 0, 0, 0, 0; -t114 * t137 + t140 * t151, -t145 * t137, 0, 0, 0, 0; 0, -t124 * t137, 0, 0, 0, 0;];
JR_rot  = t1;

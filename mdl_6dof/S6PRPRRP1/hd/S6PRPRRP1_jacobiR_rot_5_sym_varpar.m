% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRPRRP1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1,theta3]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 09:33
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6PRPRRP1_jacobiR_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP1_jacobiR_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP1_jacobiR_rot_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 09:33:29
% EndTime: 2019-02-22 09:33:29
% DurationCPUTime: 0.10s
% Computational Cost: add. (114->28), mult. (319->65), div. (0->0), fcn. (454->12), ass. (0->33)
t144 = sin(pkin(11));
t147 = cos(pkin(11));
t152 = sin(qJ(2));
t155 = cos(qJ(2));
t140 = t152 * t144 - t155 * t147;
t146 = sin(pkin(6));
t151 = sin(qJ(4));
t165 = t146 * t151;
t154 = cos(qJ(4));
t164 = t146 * t154;
t150 = sin(qJ(5));
t163 = t150 * t154;
t153 = cos(qJ(5));
t161 = t153 * t154;
t149 = cos(pkin(6));
t157 = t155 * t144 + t152 * t147;
t139 = t157 * t149;
t145 = sin(pkin(10));
t148 = cos(pkin(10));
t159 = t148 * t139 - t145 * t140;
t158 = -t145 * t139 - t148 * t140;
t156 = t140 * t149;
t138 = t157 * t146;
t137 = t140 * t146;
t135 = t138 * t154 + t149 * t151;
t134 = -t138 * t151 + t149 * t154;
t132 = t145 * t156 - t148 * t157;
t129 = -t145 * t157 - t148 * t156;
t127 = t145 * t165 + t154 * t158;
t126 = t145 * t164 - t151 * t158;
t125 = -t148 * t165 + t154 * t159;
t124 = -t148 * t164 - t151 * t159;
t1 = [0, t132 * t161 + t150 * t158, 0, t126 * t153, -t127 * t150 - t132 * t153, 0; 0, t129 * t161 + t150 * t159, 0, t124 * t153, -t125 * t150 - t129 * t153, 0; 0, -t137 * t161 + t138 * t150, 0, t134 * t153, -t135 * t150 + t137 * t153, 0; 0, -t132 * t163 + t153 * t158, 0, -t126 * t150, -t127 * t153 + t132 * t150, 0; 0, -t129 * t163 + t153 * t159, 0, -t124 * t150, -t125 * t153 + t129 * t150, 0; 0, t137 * t163 + t138 * t153, 0, -t134 * t150, -t135 * t153 - t137 * t150, 0; 0, t132 * t151, 0, t127, 0, 0; 0, t129 * t151, 0, t125, 0, 0; 0, -t137 * t151, 0, t135, 0, 0;];
JR_rot  = t1;

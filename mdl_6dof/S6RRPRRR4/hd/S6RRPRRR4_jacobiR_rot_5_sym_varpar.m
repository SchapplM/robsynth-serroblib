% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPRRR4
% Use Code from Maple symbolic Code Generation
%
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6,theta3]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 11:41
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRPRRR4_jacobiR_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR4_jacobiR_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR4_jacobiR_rot_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 11:41:30
% EndTime: 2019-02-22 11:41:30
% DurationCPUTime: 0.07s
% Computational Cost: add. (124->21), mult. (238->36), div. (0->0), fcn. (348->10), ass. (0->30)
t140 = sin(pkin(6));
t144 = sin(qJ(1));
t150 = t140 * t144;
t146 = cos(qJ(1));
t149 = t140 * t146;
t142 = cos(pkin(6));
t139 = sin(pkin(12));
t141 = cos(pkin(12));
t143 = sin(qJ(2));
t145 = cos(qJ(2));
t148 = t145 * t139 + t143 * t141;
t131 = t148 * t142;
t132 = t143 * t139 - t145 * t141;
t123 = t146 * t131 - t144 * t132;
t138 = qJ(4) + qJ(5);
t136 = sin(t138);
t137 = cos(t138);
t119 = -t123 * t137 + t136 * t149;
t125 = -t144 * t131 - t146 * t132;
t147 = t123 * t136 + t137 * t149;
t130 = t132 * t142;
t129 = t148 * t140;
t128 = t132 * t140;
t127 = -t129 * t137 - t142 * t136;
t126 = -t129 * t136 + t142 * t137;
t124 = t144 * t130 - t146 * t148;
t122 = -t146 * t130 - t144 * t148;
t121 = t125 * t137 + t136 * t150;
t120 = -t125 * t136 + t137 * t150;
t1 = [t119, t124 * t137, 0, t120, t120, 0; t121, t122 * t137, 0, -t147, -t147, 0; 0, -t128 * t137, 0, t126, t126, 0; t147, -t124 * t136, 0, -t121, -t121, 0; t120, -t122 * t136, 0, t119, t119, 0; 0, t128 * t136, 0, t127, t127, 0; t122, t125, 0, 0, 0, 0; -t124, t123, 0, 0, 0, 0; 0, t129, 0, 0, 0, 0;];
JR_rot  = t1;

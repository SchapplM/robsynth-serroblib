% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPPRR5
% Use Code from Maple symbolic Code Generation
%
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 11:15
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRPPRR5_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR5_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR5_jacobiR_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 11:15:03
% EndTime: 2019-02-22 11:15:03
% DurationCPUTime: 0.14s
% Computational Cost: add. (74->28), mult. (219->64), div. (0->0), fcn. (320->10), ass. (0->34)
t133 = sin(qJ(2));
t134 = sin(qJ(1));
t137 = cos(qJ(2));
t138 = cos(qJ(1));
t150 = cos(pkin(6));
t139 = t138 * t150;
t124 = t133 * t139 + t134 * t137;
t132 = sin(qJ(5));
t136 = cos(qJ(5));
t130 = sin(pkin(6));
t144 = t130 * t138;
t117 = t124 * t136 + t132 * t144;
t123 = t134 * t133 - t137 * t139;
t131 = sin(qJ(6));
t135 = cos(qJ(6));
t152 = t117 * t131 + t123 * t135;
t151 = -t117 * t135 + t123 * t131;
t147 = t130 * t132;
t146 = t130 * t136;
t145 = t130 * t137;
t143 = t131 * t136;
t142 = t135 * t136;
t141 = t136 * t137;
t140 = t134 * t150;
t116 = -t124 * t132 + t136 * t144;
t126 = -t133 * t140 + t138 * t137;
t125 = -t138 * t133 - t137 * t140;
t122 = -t132 * t150 + t133 * t146;
t121 = -t133 * t147 - t136 * t150;
t120 = t126 * t136 - t134 * t147;
t119 = t126 * t132 + t134 * t146;
t115 = t120 * t135 + t125 * t131;
t114 = -t120 * t131 + t125 * t135;
t1 = [t151, t125 * t142 - t126 * t131, 0, 0, -t119 * t135, t114; t115, -t123 * t142 - t124 * t131, 0, 0, t116 * t135, -t152; 0 (-t131 * t133 + t135 * t141) * t130, 0, 0, t121 * t135, -t122 * t131 + t135 * t145; t152, -t125 * t143 - t126 * t135, 0, 0, t119 * t131, -t115; t114, t123 * t143 - t124 * t135, 0, 0, -t116 * t131, t151; 0 (-t131 * t141 - t133 * t135) * t130, 0, 0, -t121 * t131, -t122 * t135 - t131 * t145; t116, t125 * t132, 0, 0, t120, 0; t119, -t123 * t132, 0, 0, t117, 0; 0, t132 * t145, 0, 0, t122, 0;];
JR_rot  = t1;

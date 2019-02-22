% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRPRPR2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3,theta5]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 09:29
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6PRPRPR2_jacobiR_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR2_jacobiR_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRPR2_jacobiR_rot_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 09:29:50
% EndTime: 2019-02-22 09:29:50
% DurationCPUTime: 0.10s
% Computational Cost: add. (74->22), mult. (211->53), div. (0->0), fcn. (300->12), ass. (0->30)
t130 = sin(pkin(12));
t140 = cos(qJ(4));
t148 = t130 * t140;
t133 = sin(pkin(6));
t138 = sin(qJ(4));
t147 = t133 * t138;
t146 = t133 * t140;
t134 = cos(pkin(12));
t145 = t134 * t140;
t137 = cos(pkin(6));
t131 = sin(pkin(11));
t135 = cos(pkin(11));
t139 = sin(qJ(2));
t141 = cos(qJ(2));
t142 = t141 * t131 + t139 * t135;
t126 = t142 * t137;
t127 = t139 * t131 - t141 * t135;
t132 = sin(pkin(10));
t136 = cos(pkin(10));
t144 = t136 * t126 - t132 * t127;
t143 = -t132 * t126 - t136 * t127;
t125 = t127 * t137;
t124 = t142 * t133;
t123 = t127 * t133;
t122 = -t124 * t138 + t137 * t140;
t120 = t132 * t125 - t136 * t142;
t117 = -t136 * t125 - t132 * t142;
t115 = t132 * t146 - t138 * t143;
t114 = -t136 * t146 - t138 * t144;
t1 = [0, t120 * t145 + t130 * t143, 0, t115 * t134, 0, 0; 0, t117 * t145 + t130 * t144, 0, t114 * t134, 0, 0; 0, -t123 * t145 + t124 * t130, 0, t122 * t134, 0, 0; 0, -t120 * t148 + t134 * t143, 0, -t115 * t130, 0, 0; 0, -t117 * t148 + t134 * t144, 0, -t114 * t130, 0, 0; 0, t123 * t148 + t124 * t134, 0, -t122 * t130, 0, 0; 0, t120 * t138, 0, t132 * t147 + t140 * t143, 0, 0; 0, t117 * t138, 0, -t136 * t147 + t140 * t144, 0, 0; 0, -t123 * t138, 0, t124 * t140 + t137 * t138, 0, 0;];
JR_rot  = t1;

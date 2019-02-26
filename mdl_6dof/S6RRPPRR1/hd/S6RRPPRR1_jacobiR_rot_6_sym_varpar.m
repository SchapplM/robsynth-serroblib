% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPPRR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta3]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:28
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRPPRR1_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR1_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR1_jacobiR_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:28:06
% EndTime: 2019-02-26 21:28:07
% DurationCPUTime: 0.09s
% Computational Cost: add. (106->21), mult. (134->24), div. (0->0), fcn. (202->8), ass. (0->26)
t123 = qJ(2) + pkin(10);
t122 = cos(t123);
t125 = sin(qJ(5));
t141 = sin(t123);
t142 = cos(qJ(5));
t117 = -t122 * t125 + t141 * t142;
t116 = t122 * t142 + t141 * t125;
t126 = sin(qJ(1));
t112 = t117 * t126;
t124 = sin(qJ(6));
t140 = t112 * t124;
t127 = cos(qJ(6));
t139 = t112 * t127;
t128 = cos(qJ(1));
t115 = t117 * t128;
t138 = t115 * t124;
t137 = t115 * t127;
t136 = t116 * t124;
t135 = t116 * t127;
t113 = t116 * t126;
t130 = -t113 * t127 - t128 * t124;
t129 = t113 * t124 - t128 * t127;
t114 = t116 * t128;
t111 = t114 * t127 - t126 * t124;
t110 = -t114 * t124 - t126 * t127;
t1 = [t130, -t137, 0, 0, t137, t110; t111, -t139, 0, 0, t139, -t129; 0, t135, 0, 0, -t135, -t117 * t124; t129, t138, 0, 0, -t138, -t111; t110, t140, 0, 0, -t140, t130; 0, -t136, 0, 0, t136, -t117 * t127; t112, -t114, 0, 0, t114, 0; -t115, -t113, 0, 0, t113, 0; 0, -t117, 0, 0, t117, 0;];
JR_rot  = t1;

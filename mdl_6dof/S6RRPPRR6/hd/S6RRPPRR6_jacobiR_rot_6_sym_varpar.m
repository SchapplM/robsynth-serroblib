% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPPRR6
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta4]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:31
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRPPRR6_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR6_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR6_jacobiR_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:31:25
% EndTime: 2019-02-26 21:31:25
% DurationCPUTime: 0.09s
% Computational Cost: add. (106->21), mult. (134->24), div. (0->0), fcn. (202->8), ass. (0->26)
t121 = pkin(10) + qJ(5);
t120 = sin(t121);
t125 = cos(qJ(2));
t139 = sin(qJ(2));
t140 = cos(t121);
t115 = -t125 * t120 + t139 * t140;
t114 = t139 * t120 + t125 * t140;
t123 = sin(qJ(1));
t110 = t115 * t123;
t122 = sin(qJ(6));
t138 = t110 * t122;
t124 = cos(qJ(6));
t137 = t110 * t124;
t126 = cos(qJ(1));
t113 = t115 * t126;
t136 = t113 * t122;
t135 = t113 * t124;
t134 = t114 * t122;
t133 = t114 * t124;
t111 = t114 * t123;
t128 = -t111 * t124 - t126 * t122;
t127 = t111 * t122 - t126 * t124;
t112 = t114 * t126;
t109 = t112 * t124 - t123 * t122;
t108 = -t112 * t122 - t123 * t124;
t1 = [t128, -t135, 0, 0, t135, t108; t109, -t137, 0, 0, t137, -t127; 0, t133, 0, 0, -t133, -t115 * t122; t127, t136, 0, 0, -t136, -t109; t108, t138, 0, 0, -t138, t128; 0, -t134, 0, 0, t134, -t115 * t124; t110, -t112, 0, 0, t112, 0; -t113, -t111, 0, 0, t111, 0; 0, -t115, 0, 0, t115, 0;];
JR_rot  = t1;

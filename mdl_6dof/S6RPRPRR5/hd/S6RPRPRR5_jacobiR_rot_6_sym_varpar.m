% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRPRR5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:51
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RPRPRR5_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR5_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR5_jacobiR_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:51:15
% EndTime: 2019-02-26 20:51:15
% DurationCPUTime: 0.08s
% Computational Cost: add. (106->21), mult. (134->24), div. (0->0), fcn. (202->8), ass. (0->26)
t120 = pkin(10) + qJ(3);
t119 = cos(t120);
t122 = sin(qJ(5));
t138 = sin(t120);
t139 = cos(qJ(5));
t114 = -t119 * t122 + t138 * t139;
t113 = t119 * t139 + t138 * t122;
t123 = sin(qJ(1));
t109 = t114 * t123;
t121 = sin(qJ(6));
t137 = t109 * t121;
t124 = cos(qJ(6));
t136 = t109 * t124;
t125 = cos(qJ(1));
t112 = t114 * t125;
t135 = t112 * t121;
t134 = t112 * t124;
t133 = t113 * t121;
t132 = t113 * t124;
t110 = t113 * t123;
t127 = -t110 * t124 - t125 * t121;
t126 = t110 * t121 - t125 * t124;
t111 = t113 * t125;
t108 = t111 * t124 - t123 * t121;
t107 = -t111 * t121 - t123 * t124;
t1 = [t127, 0, -t134, 0, t134, t107; t108, 0, -t136, 0, t136, -t126; 0, 0, t132, 0, -t132, -t114 * t121; t126, 0, t135, 0, -t135, -t108; t107, 0, t137, 0, -t137, t127; 0, 0, -t133, 0, t133, -t114 * t124; t109, 0, -t111, 0, t111, 0; -t112, 0, -t110, 0, t110, 0; 0, 0, -t114, 0, t114, 0;];
JR_rot  = t1;

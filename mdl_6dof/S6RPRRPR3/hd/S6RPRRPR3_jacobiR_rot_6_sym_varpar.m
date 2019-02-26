% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRPR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:02
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RPRRPR3_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR3_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPR3_jacobiR_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:02:02
% EndTime: 2019-02-26 21:02:02
% DurationCPUTime: 0.06s
% Computational Cost: add. (110->27), mult. (148->34), div. (0->0), fcn. (221->8), ass. (0->24)
t102 = cos(qJ(6));
t103 = cos(qJ(4));
t100 = sin(qJ(4));
t104 = cos(qJ(3));
t111 = t100 * t104;
t98 = qJ(1) + pkin(10);
t96 = sin(t98);
t97 = cos(t98);
t88 = t97 * t103 + t96 * t111;
t110 = t103 * t104;
t89 = -t97 * t100 + t96 * t110;
t99 = sin(qJ(6));
t114 = t88 * t102 - t89 * t99;
t101 = sin(qJ(3));
t107 = t100 * t99 + t102 * t103;
t113 = t107 * t101;
t108 = t89 * t102 + t88 * t99;
t90 = -t96 * t103 + t97 * t111;
t91 = t96 * t100 + t97 * t110;
t86 = t91 * t102 + t90 * t99;
t85 = t90 * t102 - t91 * t99;
t106 = t100 * t102 - t103 * t99;
t92 = t106 * t101;
t1 = [-t108, 0, -t97 * t113, -t85, 0, t85; t86, 0, -t96 * t113, -t114, 0, t114; 0, 0, t107 * t104, -t92, 0, t92; -t114, 0, -t97 * t92, t86, 0, -t86; t85, 0, -t96 * t92, t108, 0, -t108; 0, 0, t106 * t104, t113, 0, -t113; t96 * t101, 0, -t97 * t104, 0, 0, 0; -t97 * t101, 0, -t96 * t104, 0, 0, 0; 0, 0, -t101, 0, 0, 0;];
JR_rot  = t1;

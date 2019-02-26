% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRRR7
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:18
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RPRRRR7_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR7_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRR7_jacobiR_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:18:11
% EndTime: 2019-02-26 21:18:11
% DurationCPUTime: 0.04s
% Computational Cost: add. (102->23), mult. (64->20), div. (0->0), fcn. (111->6), ass. (0->24)
t104 = cos(qJ(1));
t100 = qJ(3) + qJ(4) + qJ(5);
t98 = sin(t100);
t112 = t104 * t98;
t103 = cos(qJ(6));
t111 = t98 * t103;
t101 = sin(qJ(6));
t102 = sin(qJ(1));
t110 = t102 * t101;
t109 = t102 * t103;
t108 = t104 * t101;
t107 = t104 * t103;
t99 = cos(t100);
t106 = t99 * t110;
t105 = t99 * t107;
t97 = t102 * t98;
t96 = t98 * t101;
t95 = t99 * t108;
t94 = t99 * t109;
t93 = t98 * t107 - t110;
t92 = t98 * t108 + t109;
t91 = t98 * t109 + t108;
t90 = -t98 * t110 + t107;
t1 = [t93, 0, t94, t94, t94, t90; t91, 0, -t105, -t105, -t105, t92; 0, 0, -t111, -t111, -t111, -t99 * t101; -t92, 0, -t106, -t106, -t106, -t91; t90, 0, t95, t95, t95, t93; 0, 0, t96, t96, t96, -t99 * t103; -t104 * t99, 0, t97, t97, t97, 0; -t102 * t99, 0, -t112, -t112, -t112, 0; 0, 0, t99, t99, t99, 0;];
JR_rot  = t1;

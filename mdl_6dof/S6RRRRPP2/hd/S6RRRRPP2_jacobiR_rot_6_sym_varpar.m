% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRPP2
% Use Code from Maple symbolic Code Generation
%
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:25
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRRRPP2_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP2_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRRPP2_jacobiR_rot_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:25:31
% EndTime: 2019-02-26 22:25:31
% DurationCPUTime: 0.04s
% Computational Cost: add. (54->23), mult. (52->20), div. (0->0), fcn. (90->6), ass. (0->24)
t100 = sin(qJ(1));
t98 = qJ(2) + qJ(3);
t97 = cos(t98);
t112 = t100 * t97;
t99 = sin(qJ(4));
t111 = t100 * t99;
t102 = cos(qJ(1));
t110 = t102 * t97;
t109 = t102 * t99;
t101 = cos(qJ(4));
t108 = t100 * t101;
t107 = t102 * t101;
t96 = sin(t98);
t106 = t96 * t111;
t105 = t96 * t109;
t104 = t96 * t108;
t103 = t96 * t107;
t95 = t97 * t101;
t94 = t97 * t99;
t93 = t97 * t107 + t111;
t92 = t97 * t109 - t108;
t91 = t97 * t108 - t109;
t90 = -t97 * t111 - t107;
t1 = [-t91, -t103, -t103, -t92, 0, 0; t93, -t104, -t104, t90, 0, 0; 0, t95, t95, -t96 * t99, 0, 0; t90, -t105, -t105, t93, 0, 0; t92, -t106, -t106, t91, 0, 0; 0, t94, t94, t96 * t101, 0, 0; t100 * t96, -t110, -t110, 0, 0, 0; -t102 * t96, -t112, -t112, 0, 0, 0; 0, -t96, -t96, 0, 0, 0;];
JR_rot  = t1;

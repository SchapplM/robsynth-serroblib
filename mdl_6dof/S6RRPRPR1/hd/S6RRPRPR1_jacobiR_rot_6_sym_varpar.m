% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRPR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta3,theta5]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 11:22
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRPRPR1_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR1_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR1_jacobiR_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 11:22:41
% EndTime: 2019-02-22 11:22:41
% DurationCPUTime: 0.04s
% Computational Cost: add. (107->17), mult. (52->20), div. (0->0), fcn. (90->6), ass. (0->25)
t105 = qJ(2) + pkin(10) + qJ(4);
t102 = cos(t105);
t106 = pkin(11) + qJ(6);
t103 = sin(t106);
t115 = t102 * t103;
t107 = sin(qJ(1));
t114 = t107 * t103;
t104 = cos(t106);
t113 = t107 * t104;
t108 = cos(qJ(1));
t112 = t108 * t103;
t111 = t108 * t104;
t101 = sin(t105);
t110 = t101 * t113;
t109 = t101 * t111;
t100 = t108 * t102;
t99 = t107 * t102;
t98 = t102 * t104;
t97 = t101 * t112;
t96 = t101 * t114;
t95 = t102 * t111 + t114;
t94 = -t102 * t112 + t113;
t93 = -t102 * t113 + t112;
t92 = t102 * t114 + t111;
t1 = [t93, -t109, 0, -t109, 0, t94; t95, -t110, 0, -t110, 0, -t92; 0, t98, 0, t98, 0, -t101 * t103; t92, t97, 0, t97, 0, -t95; t94, t96, 0, t96, 0, t93; 0, -t115, 0, -t115, 0, -t101 * t104; -t107 * t101, t100, 0, t100, 0, 0; t108 * t101, t99, 0, t99, 0, 0; 0, t101, 0, t101, 0, 0;];
JR_rot  = t1;

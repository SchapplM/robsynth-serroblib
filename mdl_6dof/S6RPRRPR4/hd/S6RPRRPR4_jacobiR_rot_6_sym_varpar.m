% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRPR4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2,theta5]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 10:45
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RPRRPR4_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR4_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRPR4_jacobiR_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 10:45:24
% EndTime: 2019-02-22 10:45:25
% DurationCPUTime: 0.05s
% Computational Cost: add. (107->17), mult. (52->20), div. (0->0), fcn. (90->6), ass. (0->25)
t103 = pkin(10) + qJ(3) + qJ(4);
t100 = cos(t103);
t104 = pkin(11) + qJ(6);
t101 = sin(t104);
t113 = t100 * t101;
t105 = sin(qJ(1));
t112 = t105 * t101;
t102 = cos(t104);
t111 = t105 * t102;
t106 = cos(qJ(1));
t110 = t106 * t101;
t109 = t106 * t102;
t99 = sin(t103);
t108 = t99 * t111;
t107 = t99 * t109;
t98 = t106 * t100;
t97 = t105 * t100;
t96 = t100 * t102;
t95 = t99 * t110;
t94 = t99 * t112;
t93 = t100 * t109 + t112;
t92 = -t100 * t110 + t111;
t91 = -t100 * t111 + t110;
t90 = t100 * t112 + t109;
t1 = [t91, 0, -t107, -t107, 0, t92; t93, 0, -t108, -t108, 0, -t90; 0, 0, t96, t96, 0, -t99 * t101; t90, 0, t95, t95, 0, -t93; t92, 0, t94, t94, 0, t91; 0, 0, -t113, -t113, 0, -t99 * t102; -t105 * t99, 0, t98, t98, 0, 0; t106 * t99, 0, t97, t97, 0, 0; 0, 0, t99, t99, 0, 0;];
JR_rot  = t1;

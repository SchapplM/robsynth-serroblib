% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRPR5
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
% Datum: 2019-02-22 10:46
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RPRRPR5_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR5_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPR5_jacobiR_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 10:46:02
% EndTime: 2019-02-22 10:46:02
% DurationCPUTime: 0.04s
% Computational Cost: add. (74->13), mult. (52->20), div. (0->0), fcn. (90->6), ass. (0->24)
t104 = pkin(10) + qJ(3) + qJ(4);
t102 = sin(t104);
t106 = sin(qJ(1));
t114 = t106 * t102;
t105 = sin(qJ(6));
t113 = t106 * t105;
t107 = cos(qJ(6));
t112 = t106 * t107;
t108 = cos(qJ(1));
t111 = t108 * t102;
t110 = t108 * t105;
t109 = t108 * t107;
t103 = cos(t104);
t101 = t102 * t107;
t100 = t102 * t105;
t99 = t103 * t109;
t98 = t103 * t110;
t97 = t103 * t112;
t96 = t103 * t113;
t95 = -t102 * t113 + t109;
t94 = t102 * t112 + t110;
t93 = t102 * t110 + t112;
t92 = t102 * t109 - t113;
t1 = [t95, 0, t98, t98, 0, t92; t93, 0, t96, t96, 0, t94; 0, 0, t100, t100, 0, -t103 * t107; -t94, 0, t99, t99, 0, -t93; t92, 0, t97, t97, 0, t95; 0, 0, t101, t101, 0, t103 * t105; -t106 * t103, 0, -t111, -t111, 0, 0; t108 * t103, 0, -t114, -t114, 0, 0; 0, 0, t103, t103, 0, 0;];
JR_rot  = t1;

% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRRRP1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 10:51
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RPRRRP1_jacobiR_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP1_jacobiR_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP1_jacobiR_rot_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 10:51:04
% EndTime: 2019-02-22 10:51:04
% DurationCPUTime: 0.11s
% Computational Cost: add. (77->17), mult. (52->20), div. (0->0), fcn. (90->6), ass. (0->23)
t102 = sin(qJ(5));
t101 = qJ(3) + qJ(4);
t98 = sin(t101);
t108 = t98 * t102;
t103 = cos(qJ(5));
t107 = t98 * t103;
t99 = cos(t101);
t106 = t99 * t102;
t95 = t99 * t103;
t100 = qJ(1) + pkin(10);
t96 = sin(t100);
t105 = t96 * t107;
t97 = cos(t100);
t104 = t97 * t107;
t94 = t97 * t99;
t93 = t96 * t99;
t92 = t97 * t108;
t91 = t96 * t108;
t90 = t96 * t102 + t95 * t97;
t89 = t96 * t103 - t106 * t97;
t88 = t97 * t102 - t95 * t96;
t87 = t97 * t103 + t106 * t96;
t1 = [t88, 0, -t104, -t104, t89, 0; t90, 0, -t105, -t105, -t87, 0; 0, 0, t95, t95, -t108, 0; t87, 0, t92, t92, -t90, 0; t89, 0, t91, t91, t88, 0; 0, 0, -t106, -t106, -t107, 0; -t96 * t98, 0, t94, t94, 0, 0; t97 * t98, 0, t93, t93, 0, 0; 0, 0, t98, t98, 0, 0;];
JR_rot  = t1;

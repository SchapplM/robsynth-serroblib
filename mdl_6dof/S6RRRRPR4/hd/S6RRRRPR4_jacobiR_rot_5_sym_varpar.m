% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRRPR4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6,theta5]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 12:19
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRRRPR4_jacobiR_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR4_jacobiR_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR4_jacobiR_rot_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 12:19:29
% EndTime: 2019-02-22 12:19:29
% DurationCPUTime: 0.04s
% Computational Cost: add. (77->17), mult. (52->20), div. (0->0), fcn. (90->6), ass. (0->23)
t100 = qJ(4) + pkin(11);
t96 = sin(t100);
t101 = qJ(2) + qJ(3);
t99 = cos(t101);
t108 = t99 * t96;
t102 = sin(qJ(1));
t98 = sin(t101);
t107 = t102 * t98;
t94 = t102 * t99;
t103 = cos(qJ(1));
t106 = t103 * t98;
t95 = t103 * t99;
t97 = cos(t100);
t105 = t97 * t107;
t104 = t97 * t106;
t93 = t99 * t97;
t92 = t96 * t106;
t91 = t96 * t107;
t90 = t102 * t96 + t97 * t95;
t89 = t102 * t97 - t96 * t95;
t88 = t103 * t96 - t97 * t94;
t87 = t103 * t97 + t96 * t94;
t1 = [t88, -t104, -t104, t89, 0, 0; t90, -t105, -t105, -t87, 0, 0; 0, t93, t93, -t98 * t96, 0, 0; t87, t92, t92, -t90, 0, 0; t89, t91, t91, t88, 0, 0; 0, -t108, -t108, -t98 * t97, 0, 0; -t107, t95, t95, 0, 0, 0; t106, t94, t94, 0, 0, 0; 0, t98, t98, 0, 0, 0;];
JR_rot  = t1;

% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRPRR4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6,theta4]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 12:05
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRRPRR4_jacobiR_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR4_jacobiR_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR4_jacobiR_rot_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 12:05:15
% EndTime: 2019-02-22 12:05:15
% DurationCPUTime: 0.04s
% Computational Cost: add. (77->17), mult. (52->20), div. (0->0), fcn. (90->6), ass. (0->23)
t98 = pkin(11) + qJ(5);
t94 = sin(t98);
t99 = qJ(2) + qJ(3);
t97 = cos(t99);
t106 = t97 * t94;
t100 = sin(qJ(1));
t96 = sin(t99);
t105 = t100 * t96;
t92 = t100 * t97;
t101 = cos(qJ(1));
t104 = t101 * t96;
t93 = t101 * t97;
t95 = cos(t98);
t103 = t95 * t105;
t102 = t95 * t104;
t91 = t97 * t95;
t90 = t94 * t104;
t89 = t94 * t105;
t88 = t100 * t94 + t95 * t93;
t87 = t100 * t95 - t94 * t93;
t86 = t101 * t94 - t95 * t92;
t85 = t101 * t95 + t94 * t92;
t1 = [t86, -t102, -t102, 0, t87, 0; t88, -t103, -t103, 0, -t85, 0; 0, t91, t91, 0, -t96 * t94, 0; t85, t90, t90, 0, -t88, 0; t87, t89, t89, 0, t86, 0; 0, -t106, -t106, 0, -t96 * t95, 0; -t105, t93, t93, 0, 0, 0; t104, t92, t92, 0, 0, 0; 0, t96, t96, 0, 0, 0;];
JR_rot  = t1;

% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S5RRRRR2
% Use Code from Maple symbolic Code Generation
%
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4]';
%
% Output:
% JR_rot [9x5]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-29 15:26
% Revision: 932832b1be1be80f59b7f1a581a1a8f328bdb39d (2019-03-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S5RRRRR2_jacobiR_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(2,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR2_jacobiR_rot_5_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'S5RRRRR2_jacobiR_rot_5_sym_varpar: pkin has to be [2x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-29 15:26:51
% EndTime: 2019-03-29 15:26:51
% DurationCPUTime: 0.04s
% Computational Cost: add. (98->18), mult. (66->20), div. (0->0), fcn. (114->6), ass. (0->27)
t104 = qJ(3) + qJ(4);
t100 = sin(t104);
t105 = qJ(1) + qJ(2);
t101 = sin(t105);
t115 = t101 * t100;
t106 = sin(qJ(5));
t114 = t101 * t106;
t107 = cos(qJ(5));
t113 = t101 * t107;
t102 = cos(t104);
t112 = t102 * t106;
t99 = t102 * t107;
t103 = cos(t105);
t111 = t103 * t106;
t110 = t103 * t107;
t109 = t100 * t113;
t108 = t100 * t110;
t98 = t103 * t102;
t97 = t103 * t100;
t96 = t101 * t102;
t95 = t100 * t111;
t94 = t100 * t114;
t93 = t102 * t110 + t114;
t92 = -t102 * t111 + t113;
t91 = -t101 * t99 + t111;
t90 = t101 * t112 + t110;
t1 = [t91, t91, -t108, -t108, t92; t93, t93, -t109, -t109, -t90; 0, 0, t99, t99, -t100 * t106; t90, t90, t95, t95, -t93; t92, t92, t94, t94, t91; 0, 0, -t112, -t112, -t100 * t107; -t115, -t115, t98, t98, 0; t97, t97, t96, t96, 0; 0, 0, t100, t100, 0;];
JR_rot  = t1;

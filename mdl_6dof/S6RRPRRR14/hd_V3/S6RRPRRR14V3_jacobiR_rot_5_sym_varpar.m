% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPRRR14V3
% Use Code from Maple symbolic Code Generation
%
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[dummy]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-04-12 15:12
% Revision: b693519ea345eb34ae9622239e7f1167217e9d53 (2019-04-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRPRRR14V3_jacobiR_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(1,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR14V3_jacobiR_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S6RRPRRR14V3_jacobiR_rot_5_sym_varpar: pkin has to be [1x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-04-12 15:12:04
% EndTime: 2019-04-12 15:12:04
% DurationCPUTime: 0.06s
% Computational Cost: add. (37->20), mult. (118->41), div. (0->0), fcn. (180->8), ass. (0->27)
t94 = sin(qJ(5));
t96 = sin(qJ(2));
t113 = t96 * t94;
t98 = cos(qJ(5));
t112 = t96 * t98;
t99 = cos(qJ(4));
t111 = t96 * t99;
t95 = sin(qJ(4));
t97 = sin(qJ(1));
t110 = t97 * t95;
t100 = cos(qJ(2));
t109 = t100 * t95;
t108 = t100 * t99;
t101 = cos(qJ(1));
t107 = t101 * t96;
t106 = t100 * t101;
t89 = -t101 * t95 + t97 * t108;
t105 = t97 * t112 - t89 * t94;
t104 = -t97 * t113 - t89 * t98;
t103 = t100 * t94 - t98 * t111;
t102 = t100 * t98 + t94 * t111;
t91 = t99 * t106 + t110;
t90 = t95 * t106 - t97 * t99;
t88 = -t101 * t99 - t97 * t109;
t87 = t94 * t107 + t91 * t98;
t86 = t98 * t107 - t91 * t94;
t1 = [t104, t103 * t101, 0, -t90 * t98, t86, 0; t87, t103 * t97, 0, t88 * t98, t105, 0; 0, t98 * t108 + t113, 0, -t95 * t112, -t102, 0; -t105, t102 * t101, 0, t90 * t94, -t87, 0; t86, t102 * t97, 0, -t88 * t94, t104, 0; 0, -t94 * t108 + t112, 0, t95 * t113, t103, 0; t88, -t95 * t107, 0, t91, 0, 0; t90, -t96 * t110, 0, t89, 0, 0; 0, t109, 0, t111, 0, 0;];
JR_rot  = t1;

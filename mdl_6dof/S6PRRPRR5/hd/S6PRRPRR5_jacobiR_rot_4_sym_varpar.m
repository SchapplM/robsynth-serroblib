% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6PRRPRR5
% Use Code from Maple symbolic Code Generation
%
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1,theta4]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 09:50
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6PRRPRR5_jacobiR_rot_4_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR5_jacobiR_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR5_jacobiR_rot_4_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 09:50:14
% EndTime: 2019-02-22 09:50:14
% DurationCPUTime: 0.06s
% Computational Cost: add. (34->19), mult. (107->48), div. (0->0), fcn. (156->10), ass. (0->27)
t100 = sin(qJ(3));
t96 = sin(pkin(6));
t112 = t100 * t96;
t102 = cos(qJ(3));
t94 = sin(pkin(12));
t111 = t102 * t94;
t110 = t102 * t96;
t97 = cos(pkin(12));
t109 = t102 * t97;
t101 = sin(qJ(2));
t95 = sin(pkin(11));
t108 = t95 * t101;
t103 = cos(qJ(2));
t107 = t95 * t103;
t98 = cos(pkin(11));
t106 = t98 * t101;
t105 = t98 * t103;
t104 = t102 * t103;
t99 = cos(pkin(6));
t93 = -t101 * t112 + t99 * t102;
t92 = -t99 * t108 + t105;
t91 = -t99 * t107 - t106;
t90 = t99 * t106 + t107;
t89 = t99 * t105 - t108;
t88 = -t92 * t100 + t95 * t110;
t87 = -t90 * t100 - t98 * t110;
t1 = [0, t91 * t109 + t92 * t94, t88 * t97, 0, 0, 0; 0, t89 * t109 + t90 * t94, t87 * t97, 0, 0, 0; 0 (t101 * t94 + t97 * t104) * t96, t93 * t97, 0, 0, 0; 0, -t91 * t111 + t92 * t97, -t88 * t94, 0, 0, 0; 0, -t89 * t111 + t90 * t97, -t87 * t94, 0, 0, 0; 0 (t101 * t97 - t94 * t104) * t96, -t93 * t94, 0, 0, 0; 0, t91 * t100, t92 * t102 + t95 * t112, 0, 0, 0; 0, t89 * t100, t90 * t102 - t98 * t112, 0, 0, 0; 0, t103 * t112, t99 * t100 + t101 * t110, 0, 0, 0;];
JR_rot  = t1;

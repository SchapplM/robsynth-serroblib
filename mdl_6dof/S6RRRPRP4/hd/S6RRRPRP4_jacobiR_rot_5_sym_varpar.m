% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRPRP4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 11:58
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRRPRP4_jacobiR_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP4_jacobiR_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPRP4_jacobiR_rot_5_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 11:58:02
% EndTime: 2019-02-22 11:58:02
% DurationCPUTime: 0.04s
% Computational Cost: add. (44->13), mult. (52->20), div. (0->0), fcn. (90->6), ass. (0->24)
t103 = qJ(2) + qJ(3);
t101 = sin(t103);
t105 = sin(qJ(1));
t113 = t105 * t101;
t104 = sin(qJ(5));
t112 = t105 * t104;
t106 = cos(qJ(5));
t111 = t105 * t106;
t107 = cos(qJ(1));
t110 = t107 * t101;
t109 = t107 * t104;
t108 = t107 * t106;
t102 = cos(t103);
t100 = t101 * t106;
t99 = t101 * t104;
t98 = t102 * t108;
t97 = t102 * t109;
t96 = t102 * t111;
t95 = t102 * t112;
t94 = -t101 * t112 + t108;
t93 = t101 * t111 + t109;
t92 = t101 * t109 + t111;
t91 = t101 * t108 - t112;
t1 = [t94, t97, t97, 0, t91, 0; t92, t95, t95, 0, t93, 0; 0, t99, t99, 0, -t102 * t106, 0; -t93, t98, t98, 0, -t92, 0; t91, t96, t96, 0, t94, 0; 0, t100, t100, 0, t102 * t104, 0; -t105 * t102, -t110, -t110, 0, 0, 0; t107 * t102, -t113, -t113, 0, 0, 0; 0, t102, t102, 0, 0, 0;];
JR_rot  = t1;

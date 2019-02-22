% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPRP8
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
% Datum: 2019-02-22 12:00
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRRPRP8_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP8_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPRP8_jacobiR_rot_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 12:00:30
% EndTime: 2019-02-22 12:00:30
% DurationCPUTime: 0.11s
% Computational Cost: add. (50->25), mult. (148->32), div. (0->0), fcn. (221->8), ass. (0->23)
t103 = sin(qJ(5));
t107 = cos(qJ(5));
t104 = sin(qJ(3));
t108 = cos(qJ(3));
t110 = cos(qJ(1));
t106 = sin(qJ(1));
t109 = cos(qJ(2));
t118 = t106 * t109;
t97 = t104 * t118 + t110 * t108;
t98 = -t110 * t104 + t108 * t118;
t116 = t98 * t103 - t97 * t107;
t105 = sin(qJ(2));
t112 = t103 * t104 + t107 * t108;
t120 = t112 * t105;
t117 = t110 * t109;
t115 = t97 * t103 + t98 * t107;
t100 = t106 * t104 + t108 * t117;
t99 = t104 * t117 - t106 * t108;
t93 = t100 * t107 + t99 * t103;
t114 = t100 * t103 - t99 * t107;
t113 = t103 * t108 - t104 * t107;
t95 = t113 * t105;
t1 = [-t115, -t110 * t120, t114, 0, -t114, 0; t93, -t106 * t120, t116, 0, -t116, 0; 0, t112 * t109, t95, 0, -t95, 0; t116, t110 * t95, t93, 0, -t93, 0; -t114, t106 * t95, t115, 0, -t115, 0; 0, -t113 * t109, t120, 0, -t120, 0; t106 * t105, -t117, 0, 0, 0, 0; -t110 * t105, -t118, 0, 0, 0, 0; 0, -t105, 0, 0, 0, 0;];
JR_rot  = t1;

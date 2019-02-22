% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRRPPR8
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 11:54
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRRPPR8_jacobiR_rot_4_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR8_jacobiR_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR8_jacobiR_rot_4_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 11:54:17
% EndTime: 2019-02-22 11:54:17
% DurationCPUTime: 0.05s
% Computational Cost: add. (26->14), mult. (87->31), div. (0->0), fcn. (134->8), ass. (0->25)
t102 = sin(pkin(6));
t105 = sin(qJ(2));
t119 = t102 * t105;
t107 = cos(qJ(3));
t118 = t102 * t107;
t108 = cos(qJ(2));
t117 = t102 * t108;
t109 = cos(qJ(1));
t116 = t102 * t109;
t106 = sin(qJ(1));
t115 = t106 * t105;
t114 = t106 * t108;
t113 = t109 * t105;
t112 = t109 * t108;
t104 = sin(qJ(3));
t103 = cos(pkin(6));
t99 = t103 * t113 + t114;
t111 = -t99 * t104 - t107 * t116;
t110 = t104 * t116 - t99 * t107;
t101 = -t103 * t115 + t112;
t100 = t103 * t114 + t113;
t98 = t103 * t112 - t115;
t97 = t106 * t102 * t104 + t101 * t107;
t96 = t101 * t104 - t106 * t118;
t1 = [t110, -t100 * t107, -t96, 0, 0, 0; t97, t98 * t107, t111, 0, 0, 0; 0, t107 * t117, t103 * t107 - t104 * t119, 0, 0, 0; t98, t101, 0, 0, 0, 0; t100, t99, 0, 0, 0, 0; 0, t119, 0, 0, 0, 0; t111, -t100 * t104, t97, 0, 0, 0; t96, t98 * t104, -t110, 0, 0, 0; 0, t104 * t117, t103 * t104 + t105 * t118, 0, 0, 0;];
JR_rot  = t1;

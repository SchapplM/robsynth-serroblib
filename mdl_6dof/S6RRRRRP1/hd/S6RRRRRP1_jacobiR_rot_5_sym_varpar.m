% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRRRP1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 12:26
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRRRRP1_jacobiR_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP1_jacobiR_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP1_jacobiR_rot_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 12:26:19
% EndTime: 2019-02-22 12:26:19
% DurationCPUTime: 0.09s
% Computational Cost: add. (98->19), mult. (64->20), div. (0->0), fcn. (111->6), ass. (0->24)
t106 = qJ(2) + qJ(3) + qJ(4);
t105 = cos(t106);
t107 = sin(qJ(5));
t117 = t105 * t107;
t108 = sin(qJ(1));
t116 = t108 * t107;
t109 = cos(qJ(5));
t115 = t108 * t109;
t110 = cos(qJ(1));
t114 = t110 * t107;
t113 = t110 * t109;
t104 = sin(t106);
t112 = t104 * t115;
t111 = t104 * t113;
t103 = t110 * t105;
t102 = t105 * t109;
t101 = t108 * t105;
t100 = t104 * t114;
t99 = t104 * t116;
t98 = t105 * t113 + t116;
t97 = -t105 * t114 + t115;
t96 = -t105 * t115 + t114;
t95 = t105 * t116 + t113;
t1 = [t96, -t111, -t111, -t111, t97, 0; t98, -t112, -t112, -t112, -t95, 0; 0, t102, t102, t102, -t104 * t107, 0; t95, t100, t100, t100, -t98, 0; t97, t99, t99, t99, t96, 0; 0, -t117, -t117, -t117, -t104 * t109, 0; -t108 * t104, t103, t103, t103, 0, 0; t110 * t104, t101, t101, t101, 0, 0; 0, t104, t104, t104, 0, 0;];
JR_rot  = t1;

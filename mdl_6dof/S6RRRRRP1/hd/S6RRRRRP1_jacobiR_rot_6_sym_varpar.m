% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
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

function JR_rot = S6RRRRRP1_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP1_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP1_jacobiR_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 12:26:19
% EndTime: 2019-02-22 12:26:19
% DurationCPUTime: 0.10s
% Computational Cost: add. (98->19), mult. (64->20), div. (0->0), fcn. (111->6), ass. (0->24)
t104 = qJ(2) + qJ(3) + qJ(4);
t103 = cos(t104);
t105 = sin(qJ(5));
t115 = t103 * t105;
t106 = sin(qJ(1));
t114 = t106 * t105;
t107 = cos(qJ(5));
t113 = t106 * t107;
t108 = cos(qJ(1));
t112 = t108 * t105;
t111 = t108 * t107;
t102 = sin(t104);
t110 = t102 * t113;
t109 = t102 * t111;
t101 = t108 * t103;
t100 = t103 * t107;
t99 = t106 * t103;
t98 = t102 * t112;
t97 = t102 * t114;
t96 = t103 * t111 + t114;
t95 = -t103 * t112 + t113;
t94 = -t103 * t113 + t112;
t93 = t103 * t114 + t111;
t1 = [t94, -t109, -t109, -t109, t95, 0; t96, -t110, -t110, -t110, -t93, 0; 0, t100, t100, t100, -t102 * t105, 0; t93, t98, t98, t98, -t96, 0; t95, t97, t97, t97, t94, 0; 0, -t115, -t115, -t115, -t102 * t107, 0; -t106 * t102, t101, t101, t101, 0, 0; t108 * t102, t99, t99, t99, 0, 0; 0, t102, t102, t102, 0, 0;];
JR_rot  = t1;

% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRRP8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 10:55
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RPRRRP8_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP8_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRRP8_jacobiR_rot_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 10:55:12
% EndTime: 2019-02-22 10:55:12
% DurationCPUTime: 0.04s
% Computational Cost: add. (51->20), mult. (52->20), div. (0->0), fcn. (90->6), ass. (0->24)
t105 = qJ(3) + qJ(4);
t103 = sin(t105);
t106 = sin(qJ(5));
t118 = t103 * t106;
t108 = cos(qJ(5));
t117 = t103 * t108;
t107 = sin(qJ(1));
t116 = t107 * t106;
t115 = t107 * t108;
t109 = cos(qJ(1));
t114 = t109 * t103;
t113 = t109 * t106;
t112 = t109 * t108;
t104 = cos(t105);
t111 = t104 * t113;
t110 = t104 * t112;
t102 = t107 * t103;
t101 = t104 * t115;
t100 = t104 * t116;
t99 = t103 * t112 - t116;
t98 = t103 * t113 + t115;
t97 = t103 * t115 + t113;
t96 = t103 * t116 - t112;
t1 = [t99, 0, t101, t101, -t96, 0; t97, 0, -t110, -t110, t98, 0; 0, 0, -t117, -t117, -t104 * t106, 0; -t109 * t104, 0, t102, t102, 0, 0; -t107 * t104, 0, -t114, -t114, 0, 0; 0, 0, t104, t104, 0, 0; t98, 0, t100, t100, t97, 0; t96, 0, -t111, -t111, -t99, 0; 0, 0, -t118, -t118, t104 * t108, 0;];
JR_rot  = t1;

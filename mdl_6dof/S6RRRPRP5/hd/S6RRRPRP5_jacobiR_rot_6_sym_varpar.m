% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPRP5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,theta4]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 11:58
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRRPRP5_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP5_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP5_jacobiR_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 11:58:42
% EndTime: 2019-02-22 11:58:42
% DurationCPUTime: 0.07s
% Computational Cost: add. (87->15), mult. (54->20), div. (0->0), fcn. (93->6), ass. (0->19)
t107 = qJ(3) + pkin(10) + qJ(5);
t105 = sin(t107);
t108 = sin(qJ(2));
t117 = t108 * t105;
t109 = sin(qJ(1));
t116 = t109 * t108;
t110 = cos(qJ(2));
t115 = t110 * t105;
t106 = cos(t107);
t114 = t110 * t106;
t111 = cos(qJ(1));
t113 = t111 * t108;
t112 = t111 * t110;
t103 = t108 * t106;
t102 = t109 * t105 + t106 * t112;
t101 = t105 * t112 - t109 * t106;
t100 = -t111 * t105 + t109 * t114;
t99 = -t111 * t106 - t109 * t115;
t1 = [-t100, -t106 * t113, -t101, 0, -t101, 0; t102, -t106 * t116, t99, 0, t99, 0; 0, t114, -t117, 0, -t117, 0; -t116, t112, 0, 0, 0, 0; t113, t109 * t110, 0, 0, 0, 0; 0, t108, 0, 0, 0, 0; t99, -t105 * t113, t102, 0, t102, 0; t101, -t105 * t116, t100, 0, t100, 0; 0, t115, t103, 0, t103, 0;];
JR_rot  = t1;

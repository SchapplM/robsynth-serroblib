% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPPRR8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta3]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 11:17
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRPPRR8_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR8_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR8_jacobiR_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 11:17:02
% EndTime: 2019-02-22 11:17:02
% DurationCPUTime: 0.12s
% Computational Cost: add. (94->26), mult. (148->32), div. (0->0), fcn. (221->8), ass. (0->24)
t106 = sin(qJ(2));
t103 = qJ(5) + qJ(6);
t101 = sin(t103);
t102 = cos(t103);
t104 = sin(pkin(10));
t105 = cos(pkin(10));
t112 = t101 * t105 - t102 * t104;
t116 = t106 * t112;
t107 = sin(qJ(1));
t108 = cos(qJ(2));
t115 = t107 * t108;
t109 = cos(qJ(1));
t114 = t109 * t108;
t96 = t104 * t115 + t109 * t105;
t97 = -t109 * t104 + t105 * t115;
t113 = t97 * t101 - t96 * t102;
t91 = -t96 * t101 - t97 * t102;
t111 = t101 * t104 + t102 * t105;
t95 = t111 * t106;
t99 = t107 * t104 + t105 * t114;
t98 = t104 * t114 - t107 * t105;
t93 = t98 * t101 + t99 * t102;
t92 = -t99 * t101 + t98 * t102;
t1 = [t91, -t109 * t95, 0, 0, t92, t92; t93, -t107 * t95, 0, 0, -t113, -t113; 0, t111 * t108, 0, 0, -t116, -t116; t113, t109 * t116, 0, 0, -t93, -t93; t92, t107 * t116, 0, 0, t91, t91; 0, -t112 * t108, 0, 0, -t95, -t95; t107 * t106, -t114, 0, 0, 0, 0; -t109 * t106, -t115, 0, 0, 0, 0; 0, -t106, 0, 0, 0, 0;];
JR_rot  = t1;

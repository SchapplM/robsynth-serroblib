% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRPR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta3]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:38
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRPRPR2_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR2_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR2_jacobiR_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:38:16
% EndTime: 2019-02-26 21:38:16
% DurationCPUTime: 0.04s
% Computational Cost: add. (74->13), mult. (52->20), div. (0->0), fcn. (90->6), ass. (0->24)
t107 = qJ(2) + pkin(10) + qJ(4);
t105 = sin(t107);
t109 = sin(qJ(1));
t117 = t109 * t105;
t108 = sin(qJ(6));
t116 = t109 * t108;
t110 = cos(qJ(6));
t115 = t109 * t110;
t111 = cos(qJ(1));
t114 = t111 * t105;
t113 = t111 * t108;
t112 = t111 * t110;
t106 = cos(t107);
t104 = t105 * t110;
t103 = t105 * t108;
t102 = t106 * t112;
t101 = t106 * t113;
t100 = t106 * t115;
t99 = t106 * t116;
t98 = -t105 * t116 + t112;
t97 = t105 * t115 + t113;
t96 = t105 * t113 + t115;
t95 = t105 * t112 - t116;
t1 = [t98, t101, 0, t101, 0, t95; t96, t99, 0, t99, 0, t97; 0, t103, 0, t103, 0, -t106 * t110; -t97, t102, 0, t102, 0, -t96; t95, t100, 0, t100, 0, t98; 0, t104, 0, t104, 0, t106 * t108; -t109 * t106, -t114, 0, -t114, 0, 0; t111 * t106, -t117, 0, -t117, 0, 0; 0, t106, 0, t106, 0, 0;];
JR_rot  = t1;

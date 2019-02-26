% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRRP3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:09
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RPRRRP3_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP3_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP3_jacobiR_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:09:03
% EndTime: 2019-02-26 21:09:03
% DurationCPUTime: 0.04s
% Computational Cost: add. (85->16), mult. (54->20), div. (0->0), fcn. (93->6), ass. (0->17)
t110 = qJ(4) + qJ(5);
t107 = sin(t110);
t111 = sin(qJ(3));
t115 = t111 * t107;
t108 = cos(t110);
t104 = t111 * t108;
t112 = cos(qJ(3));
t114 = t112 * t107;
t113 = t112 * t108;
t109 = qJ(1) + pkin(10);
t106 = cos(t109);
t105 = sin(t109);
t102 = t105 * t107 + t106 * t113;
t101 = -t105 * t108 + t106 * t114;
t100 = t105 * t113 - t106 * t107;
t99 = -t105 * t114 - t106 * t108;
t1 = [-t100, 0, -t106 * t104, -t101, -t101, 0; t102, 0, -t105 * t104, t99, t99, 0; 0, 0, t113, -t115, -t115, 0; -t105 * t111, 0, t106 * t112, 0, 0, 0; t106 * t111, 0, t105 * t112, 0, 0, 0; 0, 0, t111, 0, 0, 0; t99, 0, -t106 * t115, t102, t102, 0; t101, 0, -t105 * t115, t100, t100, 0; 0, 0, t114, t104, t104, 0;];
JR_rot  = t1;

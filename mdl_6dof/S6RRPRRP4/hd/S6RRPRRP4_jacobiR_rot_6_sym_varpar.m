% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRRP4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,theta3]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 11:32
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRPRRP4_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP4_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP4_jacobiR_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 11:32:44
% EndTime: 2019-02-22 11:32:44
% DurationCPUTime: 0.04s
% Computational Cost: add. (80->16), mult. (54->20), div. (0->0), fcn. (93->6), ass. (0->19)
t118 = qJ(2) + pkin(10);
t114 = sin(t118);
t119 = qJ(4) + qJ(5);
t116 = sin(t119);
t126 = t114 * t116;
t120 = sin(qJ(1));
t125 = t120 * t116;
t117 = cos(t119);
t124 = t120 * t117;
t121 = cos(qJ(1));
t123 = t121 * t116;
t122 = t121 * t117;
t115 = cos(t118);
t112 = t114 * t117;
t111 = t115 * t122 + t125;
t110 = t115 * t123 - t124;
t109 = t115 * t124 - t123;
t108 = -t115 * t125 - t122;
t1 = [-t109, -t114 * t122, 0, -t110, -t110, 0; t111, -t114 * t124, 0, t108, t108, 0; 0, t115 * t117, 0, -t126, -t126, 0; -t120 * t114, t121 * t115, 0, 0, 0, 0; t121 * t114, t120 * t115, 0, 0, 0, 0; 0, t114, 0, 0, 0, 0; t108, -t114 * t123, 0, t111, t111, 0; t110, -t114 * t125, 0, t109, t109, 0; 0, t115 * t116, 0, t112, t112, 0;];
JR_rot  = t1;

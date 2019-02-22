% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRRP6
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
% Datum: 2019-02-22 12:29
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRRRRP6_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP6_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP6_jacobiR_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 12:29:12
% EndTime: 2019-02-22 12:29:12
% DurationCPUTime: 0.04s
% Computational Cost: add. (112->17), mult. (68->20), div. (0->0), fcn. (117->6), ass. (0->19)
t111 = qJ(3) + qJ(4) + qJ(5);
t109 = sin(t111);
t112 = sin(qJ(2));
t121 = t112 * t109;
t113 = sin(qJ(1));
t120 = t113 * t112;
t114 = cos(qJ(2));
t119 = t114 * t109;
t110 = cos(t111);
t118 = t114 * t110;
t115 = cos(qJ(1));
t117 = t115 * t112;
t116 = t115 * t114;
t107 = t112 * t110;
t106 = t113 * t109 + t110 * t116;
t105 = t109 * t116 - t113 * t110;
t104 = -t115 * t109 + t113 * t118;
t103 = -t115 * t110 - t113 * t119;
t1 = [-t104, -t110 * t117, -t105, -t105, -t105, 0; t106, -t110 * t120, t103, t103, t103, 0; 0, t118, -t121, -t121, -t121, 0; -t120, t116, 0, 0, 0, 0; t117, t113 * t114, 0, 0, 0, 0; 0, t112, 0, 0, 0, 0; t103, -t109 * t117, t106, t106, t106, 0; t105, -t109 * t120, t104, t104, t104, 0; 0, t119, t107, t107, t107, 0;];
JR_rot  = t1;

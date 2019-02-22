% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRRPP6
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 12:15
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRRRPP6_jacobiR_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP6_jacobiR_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRRPP6_jacobiR_rot_5_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 12:15:12
% EndTime: 2019-02-22 12:15:12
% DurationCPUTime: 0.04s
% Computational Cost: add. (50->11), mult. (54->20), div. (0->0), fcn. (93->6), ass. (0->19)
t104 = sin(qJ(2));
t105 = sin(qJ(1));
t112 = t105 * t104;
t103 = qJ(3) + qJ(4);
t101 = sin(t103);
t106 = cos(qJ(2));
t111 = t106 * t101;
t102 = cos(t103);
t110 = t106 * t102;
t107 = cos(qJ(1));
t109 = t107 * t104;
t108 = t107 * t106;
t100 = t104 * t102;
t99 = t104 * t101;
t98 = t105 * t101 + t102 * t108;
t97 = t101 * t108 - t105 * t102;
t96 = -t107 * t101 + t105 * t110;
t95 = t107 * t102 + t105 * t111;
t1 = [-t112, t108, 0, 0, 0, 0; t109, t105 * t106, 0, 0, 0, 0; 0, t104, 0, 0, 0, 0; t96, t102 * t109, t97, t97, 0, 0; -t98, t102 * t112, t95, t95, 0, 0; 0, -t110, t99, t99, 0, 0; -t95, -t101 * t109, t98, t98, 0, 0; t97, -t101 * t112, t96, t96, 0, 0; 0, t111, t100, t100, 0, 0;];
JR_rot  = t1;

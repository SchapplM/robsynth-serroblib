% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S5RRPRR1
% Use Code from Maple symbolic Code Generation
%
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,d4,d5]';
%
% Output:
% JR_rot [9x5]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-07-18 17:22
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S5RRPRR1_jacobiR_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(4,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR1_jacobiR_rot_5_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S5RRPRR1_jacobiR_rot_5_sym_varpar: pkin has to be [4x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 17:22:52
% EndTime: 2019-07-18 17:22:52
% DurationCPUTime: 0.04s
% Computational Cost: add. (47->16), mult. (52->20), div. (0->0), fcn. (90->6), ass. (0->24)
t100 = sin(qJ(5));
t99 = qJ(2) + qJ(4);
t98 = cos(t99);
t110 = t98 * t100;
t101 = sin(qJ(1));
t109 = t101 * t100;
t102 = cos(qJ(5));
t108 = t101 * t102;
t103 = cos(qJ(1));
t107 = t103 * t100;
t106 = t103 * t102;
t97 = sin(t99);
t105 = t97 * t108;
t104 = t97 * t106;
t96 = t103 * t98;
t95 = t98 * t102;
t94 = t101 * t98;
t93 = t97 * t107;
t92 = t97 * t109;
t91 = t98 * t106 + t109;
t90 = -t98 * t107 + t108;
t89 = -t98 * t108 + t107;
t88 = t98 * t109 + t106;
t1 = [t89, -t104, 0, -t104, t90; t91, -t105, 0, -t105, -t88; 0, t95, 0, t95, -t97 * t100; t88, t93, 0, t93, -t91; t90, t92, 0, t92, t89; 0, -t110, 0, -t110, -t97 * t102; -t101 * t97, t96, 0, t96, 0; t103 * t97, t94, 0, t94, 0; 0, t97, 0, t97, 0;];
JR_rot  = t1;

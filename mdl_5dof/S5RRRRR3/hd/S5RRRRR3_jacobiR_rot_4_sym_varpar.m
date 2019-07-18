% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S5RRRRR3
% Use Code from Maple symbolic Code Generation
%
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,a5,d1,d4]';
%
% Output:
% JR_rot [9x5]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-07-18 17:19
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S5RRRRR3_jacobiR_rot_4_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR3_jacobiR_rot_4_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S5RRRRR3_jacobiR_rot_4_sym_varpar: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 17:19:38
% EndTime: 2019-07-18 17:19:38
% DurationCPUTime: 0.08s
% Computational Cost: add. (47->16), mult. (52->20), div. (0->0), fcn. (90->6), ass. (0->24)
t100 = sin(qJ(4));
t99 = qJ(2) + qJ(3);
t98 = cos(t99);
t110 = t98 * t100;
t101 = sin(qJ(1));
t109 = t101 * t100;
t102 = cos(qJ(4));
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
t91 = t106 * t98 + t109;
t90 = -t107 * t98 + t108;
t89 = -t108 * t98 + t107;
t88 = t109 * t98 + t106;
t1 = [t89, -t104, -t104, t90, 0; t91, -t105, -t105, -t88, 0; 0, t95, t95, -t97 * t100, 0; t88, t93, t93, -t91, 0; t90, t92, t92, t89, 0; 0, -t110, -t110, -t97 * t102, 0; -t101 * t97, t96, t96, 0, 0; t103 * t97, t94, t94, 0, 0; 0, t97, t97, 0, 0;];
JR_rot  = t1;

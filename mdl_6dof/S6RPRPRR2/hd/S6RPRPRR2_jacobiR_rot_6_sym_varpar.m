% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRPRR2
% Use Code from Maple symbolic Code Generation
%
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2,theta4]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:49
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RPRPRR2_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR2_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPRR2_jacobiR_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:49:39
% EndTime: 2019-02-26 20:49:39
% DurationCPUTime: 0.04s
% Computational Cost: add. (113->19), mult. (54->20), div. (0->0), fcn. (93->6), ass. (0->22)
t100 = qJ(3) + pkin(11);
t94 = sin(t100);
t102 = qJ(5) + qJ(6);
t98 = sin(t102);
t110 = t94 * t98;
t99 = cos(t102);
t109 = t94 * t99;
t101 = qJ(1) + pkin(10);
t95 = sin(t101);
t108 = t95 * t98;
t107 = t95 * t99;
t96 = cos(t100);
t106 = t96 * t98;
t105 = t96 * t99;
t97 = cos(t101);
t104 = t97 * t98;
t103 = t97 * t99;
t93 = t96 * t103 + t108;
t92 = -t96 * t104 + t107;
t91 = -t95 * t105 + t104;
t90 = t95 * t106 + t103;
t1 = [t91, 0, -t94 * t103, 0, t92, t92; t93, 0, -t94 * t107, 0, -t90, -t90; 0, 0, t105, 0, -t110, -t110; t90, 0, t94 * t104, 0, -t93, -t93; t92, 0, t94 * t108, 0, t91, t91; 0, 0, -t106, 0, -t109, -t109; -t95 * t94, 0, t97 * t96, 0, 0, 0; t97 * t94, 0, t95 * t96, 0, 0, 0; 0, 0, t94, 0, 0, 0;];
JR_rot  = t1;

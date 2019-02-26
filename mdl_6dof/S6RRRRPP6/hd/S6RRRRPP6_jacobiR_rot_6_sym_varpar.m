% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
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
% Datum: 2019-02-26 22:28
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRRRPP6_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP6_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRRPP6_jacobiR_rot_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:28:14
% EndTime: 2019-02-26 22:28:14
% DurationCPUTime: 0.04s
% Computational Cost: add. (53->15), mult. (54->20), div. (0->0), fcn. (93->6), ass. (0->19)
t102 = sin(qJ(2));
t101 = qJ(3) + qJ(4);
t99 = sin(t101);
t111 = t102 * t99;
t103 = sin(qJ(1));
t110 = t103 * t102;
t104 = cos(qJ(2));
t109 = t103 * t104;
t100 = cos(t101);
t108 = t104 * t100;
t105 = cos(qJ(1));
t107 = t105 * t102;
t106 = t105 * t104;
t97 = t102 * t100;
t96 = t100 * t106 + t103 * t99;
t95 = -t103 * t100 + t99 * t106;
t94 = t103 * t108 - t105 * t99;
t93 = -t105 * t100 - t99 * t109;
t1 = [-t110, t106, 0, 0, 0, 0; t107, t109, 0, 0, 0, 0; 0, t102, 0, 0, 0, 0; t93, -t99 * t107, t96, t96, 0, 0; t95, -t99 * t110, t94, t94, 0, 0; 0, t104 * t99, t97, t97, 0, 0; -t94, -t100 * t107, -t95, -t95, 0, 0; t96, -t100 * t110, t93, t93, 0, 0; 0, t108, -t111, -t111, 0, 0;];
JR_rot  = t1;

% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRPP5
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
% Datum: 2019-02-26 22:27
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRRRPP5_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP5_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRRPP5_jacobiR_rot_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:27:38
% EndTime: 2019-02-26 22:27:39
% DurationCPUTime: 0.04s
% Computational Cost: add. (56->18), mult. (54->20), div. (0->0), fcn. (93->6), ass. (0->19)
t96 = qJ(3) + qJ(4);
t94 = sin(t96);
t97 = sin(qJ(2));
t106 = t97 * t94;
t98 = sin(qJ(1));
t105 = t98 * t97;
t99 = cos(qJ(2));
t104 = t99 * t94;
t95 = cos(t96);
t103 = t99 * t95;
t100 = cos(qJ(1));
t102 = t100 * t97;
t101 = t100 * t99;
t92 = t97 * t95;
t91 = t95 * t101 + t98 * t94;
t90 = t94 * t101 - t98 * t95;
t89 = -t100 * t94 + t98 * t103;
t88 = -t100 * t95 - t98 * t104;
t1 = [-t89, -t95 * t102, -t90, -t90, 0, 0; t91, -t95 * t105, t88, t88, 0, 0; 0, t103, -t106, -t106, 0, 0; t88, -t94 * t102, t91, t91, 0, 0; t90, -t94 * t105, t89, t89, 0, 0; 0, t104, t92, t92, 0, 0; t105, -t101, 0, 0, 0, 0; -t102, -t98 * t99, 0, 0, 0, 0; 0, -t97, 0, 0, 0, 0;];
JR_rot  = t1;

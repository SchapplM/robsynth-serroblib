% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S5PRRRR1
% Use Code from Maple symbolic Code Generation
%
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4]';
%
% Output:
% JR_rot [9x5]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-07-18 13:29
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S5PRRRR1_jacobiR_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(2,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR1_jacobiR_rot_5_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'S5PRRRR1_jacobiR_rot_5_sym_varpar: pkin has to be [2x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 13:29:14
% EndTime: 2019-07-18 13:29:14
% DurationCPUTime: 0.08s
% Computational Cost: add. (47->16), mult. (52->20), div. (0->0), fcn. (90->6), ass. (0->24)
t93 = qJ(3) + qJ(4);
t92 = cos(t93);
t96 = cos(qJ(5));
t104 = t92 * t96;
t94 = sin(qJ(5));
t95 = sin(qJ(2));
t103 = t95 * t94;
t102 = t95 * t96;
t97 = cos(qJ(2));
t101 = t97 * t94;
t100 = t97 * t96;
t91 = sin(t93);
t99 = t91 * t102;
t98 = t91 * t100;
t90 = t97 * t92;
t89 = t92 * t94;
t88 = t95 * t92;
t87 = t91 * t101;
t86 = t91 * t103;
t85 = t100 * t92 + t103;
t84 = -t101 * t92 + t102;
t83 = -t102 * t92 + t101;
t82 = t103 * t92 + t100;
t1 = [0, t83, -t98, -t98, t84; 0, 0, -t104, -t104, t91 * t94; 0, t85, -t99, -t99, -t82; 0, t82, t87, t87, -t85; 0, 0, t89, t89, t91 * t96; 0, t84, t86, t86, t83; 0, -t95 * t91, t90, t90, 0; 0, 0, -t91, -t91, 0; 0, t97 * t91, t88, t88, 0;];
JR_rot  = t1;

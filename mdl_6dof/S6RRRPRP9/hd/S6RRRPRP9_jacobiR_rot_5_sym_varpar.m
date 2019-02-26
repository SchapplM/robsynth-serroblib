% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRPRP9
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:13
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRRPRP9_jacobiR_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP9_jacobiR_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPRP9_jacobiR_rot_5_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:13:43
% EndTime: 2019-02-26 22:13:43
% DurationCPUTime: 0.06s
% Computational Cost: add. (50->24), mult. (148->32), div. (0->0), fcn. (221->8), ass. (0->23)
t93 = sin(qJ(5));
t94 = sin(qJ(3));
t97 = cos(qJ(5));
t98 = cos(qJ(3));
t102 = t93 * t94 + t97 * t98;
t95 = sin(qJ(2));
t109 = t102 * t95;
t100 = cos(qJ(1));
t96 = sin(qJ(1));
t99 = cos(qJ(2));
t107 = t96 * t99;
t87 = t100 * t98 + t94 * t107;
t88 = -t100 * t94 + t98 * t107;
t105 = -t87 * t97 + t88 * t93;
t106 = t100 * t99;
t104 = t87 * t93 + t88 * t97;
t89 = t94 * t106 - t96 * t98;
t90 = t98 * t106 + t96 * t94;
t82 = t89 * t97 - t90 * t93;
t83 = t89 * t93 + t90 * t97;
t103 = t93 * t98 - t94 * t97;
t85 = t103 * t95;
t1 = [-t104, -t100 * t109, -t82, 0, t82, 0; t83, -t96 * t109, t105, 0, -t105, 0; 0, t102 * t99, t85, 0, -t85, 0; t105, t100 * t85, t83, 0, -t83, 0; t82, t96 * t85, t104, 0, -t104, 0; 0, -t103 * t99, t109, 0, -t109, 0; t96 * t95, -t106, 0, 0, 0, 0; -t100 * t95, -t107, 0, 0, 0, 0; 0, -t95, 0, 0, 0, 0;];
JR_rot  = t1;

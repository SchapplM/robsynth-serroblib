% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRPPR4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta4]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:41
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RPRPPR4_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR4_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR4_jacobiR_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:40:54
% EndTime: 2019-02-26 20:40:54
% DurationCPUTime: 0.05s
% Computational Cost: add. (73->22), mult. (108->32), div. (0->0), fcn. (161->8), ass. (0->26)
t87 = pkin(9) + qJ(3);
t85 = sin(t87);
t88 = sin(pkin(10));
t89 = cos(pkin(10));
t90 = sin(qJ(6));
t92 = cos(qJ(6));
t97 = t88 * t92 - t89 * t90;
t95 = t97 * t85;
t91 = sin(qJ(1));
t103 = t91 * t88;
t102 = t91 * t89;
t93 = cos(qJ(1));
t101 = t93 * t88;
t100 = t93 * t89;
t86 = cos(t87);
t80 = t86 * t103 + t100;
t81 = t86 * t102 - t101;
t99 = t80 * t92 - t81 * t90;
t98 = -t80 * t90 - t81 * t92;
t96 = t88 * t90 + t89 * t92;
t94 = t96 * t85;
t83 = t86 * t100 + t103;
t82 = t86 * t101 - t102;
t79 = t82 * t90 + t83 * t92;
t78 = t82 * t92 - t83 * t90;
t1 = [t98, 0, -t93 * t94, 0, 0, t78; t79, 0, -t91 * t94, 0, 0, t99; 0, 0, t96 * t86, 0, 0, t95; -t99, 0, -t93 * t95, 0, 0, -t79; t78, 0, -t91 * t95, 0, 0, t98; 0, 0, t97 * t86, 0, 0, -t94; t91 * t85, 0, -t93 * t86, 0, 0, 0; -t93 * t85, 0, -t91 * t86, 0, 0, 0; 0, 0, -t85, 0, 0, 0;];
JR_rot  = t1;

% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
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
% Datum: 2019-02-26 22:42
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRRRRP6_jacobiR_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP6_jacobiR_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP6_jacobiR_rot_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:42:35
% EndTime: 2019-02-26 22:42:35
% DurationCPUTime: 0.04s
% Computational Cost: add. (114->21), mult. (68->20), div. (0->0), fcn. (117->6), ass. (0->19)
t87 = qJ(3) + qJ(4) + qJ(5);
t85 = sin(t87);
t88 = sin(qJ(2));
t98 = t88 * t85;
t86 = cos(t87);
t97 = t88 * t86;
t89 = sin(qJ(1));
t96 = t89 * t88;
t90 = cos(qJ(2));
t95 = t90 * t85;
t94 = t90 * t86;
t91 = cos(qJ(1));
t93 = t91 * t88;
t92 = t91 * t90;
t84 = t89 * t85 + t86 * t92;
t83 = -t85 * t92 + t89 * t86;
t82 = t91 * t85 - t89 * t94;
t81 = t91 * t86 + t89 * t95;
t1 = [t82, -t86 * t93, t83, t83, t83, 0; t84, -t86 * t96, -t81, -t81, -t81, 0; 0, t94, -t98, -t98, -t98, 0; t81, t85 * t93, -t84, -t84, -t84, 0; t83, t85 * t96, t82, t82, t82, 0; 0, -t95, -t97, -t97, -t97, 0; -t96, t92, 0, 0, 0, 0; t93, t89 * t90, 0, 0, 0, 0; 0, t88, 0, 0, 0, 0;];
JR_rot  = t1;

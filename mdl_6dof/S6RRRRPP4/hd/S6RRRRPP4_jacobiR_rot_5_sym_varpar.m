% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRRPP4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,theta5]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:27
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRRRPP4_jacobiR_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP4_jacobiR_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPP4_jacobiR_rot_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:26:59
% EndTime: 2019-02-26 22:26:59
% DurationCPUTime: 0.04s
% Computational Cost: add. (88->17), mult. (54->20), div. (0->0), fcn. (93->6), ass. (0->19)
t85 = qJ(3) + qJ(4) + pkin(10);
t83 = sin(t85);
t86 = sin(qJ(2));
t96 = t86 * t83;
t84 = cos(t85);
t95 = t86 * t84;
t87 = sin(qJ(1));
t94 = t87 * t86;
t88 = cos(qJ(2));
t93 = t88 * t83;
t92 = t88 * t84;
t89 = cos(qJ(1));
t91 = t89 * t86;
t90 = t89 * t88;
t82 = t87 * t83 + t84 * t90;
t81 = -t83 * t90 + t87 * t84;
t80 = t89 * t83 - t87 * t92;
t79 = t89 * t84 + t87 * t93;
t1 = [t80, -t84 * t91, t81, t81, 0, 0; t82, -t84 * t94, -t79, -t79, 0, 0; 0, t92, -t96, -t96, 0, 0; t79, t83 * t91, -t82, -t82, 0, 0; t81, t83 * t94, t80, t80, 0, 0; 0, -t93, -t95, -t95, 0, 0; -t94, t90, 0, 0, 0, 0; t91, t87 * t88, 0, 0, 0, 0; 0, t86, 0, 0, 0, 0;];
JR_rot  = t1;

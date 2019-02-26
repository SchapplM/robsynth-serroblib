% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRRR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:16
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RPRRRR3_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR3_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR3_jacobiR_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:16:00
% EndTime: 2019-02-26 21:16:00
% DurationCPUTime: 0.04s
% Computational Cost: add. (154->22), mult. (68->20), div. (0->0), fcn. (117->6), ass. (0->17)
t91 = qJ(4) + qJ(5) + qJ(6);
t87 = sin(t91);
t93 = sin(qJ(3));
t98 = t93 * t87;
t88 = cos(t91);
t97 = t93 * t88;
t94 = cos(qJ(3));
t96 = t94 * t87;
t95 = t94 * t88;
t92 = qJ(1) + pkin(11);
t90 = cos(t92);
t89 = sin(t92);
t86 = t89 * t87 + t90 * t95;
t85 = t89 * t88 - t90 * t96;
t84 = t90 * t87 - t89 * t95;
t83 = t90 * t88 + t89 * t96;
t1 = [t84, 0, -t90 * t97, t85, t85, t85; t86, 0, -t89 * t97, -t83, -t83, -t83; 0, 0, t95, -t98, -t98, -t98; t83, 0, t90 * t98, -t86, -t86, -t86; t85, 0, t89 * t98, t84, t84, t84; 0, 0, -t96, -t97, -t97, -t97; -t89 * t93, 0, t90 * t94, 0, 0, 0; t90 * t93, 0, t89 * t94, 0, 0, 0; 0, 0, t93, 0, 0, 0;];
JR_rot  = t1;

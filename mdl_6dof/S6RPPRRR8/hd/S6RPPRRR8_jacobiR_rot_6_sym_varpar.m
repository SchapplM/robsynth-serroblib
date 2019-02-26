% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPPRRR8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta3]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:38
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RPPRRR8_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR8_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR8_jacobiR_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:38:39
% EndTime: 2019-02-26 20:38:39
% DurationCPUTime: 0.04s
% Computational Cost: add. (83->19), mult. (54->20), div. (0->0), fcn. (93->6), ass. (0->19)
t91 = pkin(10) + qJ(4);
t88 = cos(t91);
t92 = qJ(5) + qJ(6);
t89 = sin(t92);
t100 = t88 * t89;
t90 = cos(t92);
t99 = t88 * t90;
t93 = sin(qJ(1));
t98 = t93 * t89;
t97 = t93 * t90;
t94 = cos(qJ(1));
t96 = t94 * t89;
t95 = t94 * t90;
t87 = sin(t91);
t86 = t87 * t95 - t98;
t85 = t87 * t96 + t97;
t84 = t87 * t97 + t96;
t83 = -t87 * t98 + t95;
t1 = [t86, 0, 0, t88 * t97, t83, t83; t84, 0, 0, -t88 * t95, t85, t85; 0, 0, 0, -t87 * t90, -t100, -t100; -t85, 0, 0, -t88 * t98, -t84, -t84; t83, 0, 0, t88 * t96, t86, t86; 0, 0, 0, t87 * t89, -t99, -t99; -t94 * t88, 0, 0, t93 * t87, 0, 0; -t93 * t88, 0, 0, -t94 * t87, 0, 0; 0, 0, 0, t88, 0, 0;];
JR_rot  = t1;

% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRRPRP7
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,theta4]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:12
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRRPRP7_jacobiR_rot_4_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP7_jacobiR_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRP7_jacobiR_rot_4_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:12:38
% EndTime: 2019-02-26 22:12:38
% DurationCPUTime: 0.04s
% Computational Cost: add. (55->16), mult. (87->30), div. (0->0), fcn. (134->8), ass. (0->26)
t91 = sin(pkin(6));
t93 = sin(qJ(2));
t106 = t91 * t93;
t94 = sin(qJ(1));
t105 = t91 * t94;
t95 = cos(qJ(2));
t104 = t91 * t95;
t96 = cos(qJ(1));
t103 = t91 * t96;
t102 = t94 * t93;
t101 = t94 * t95;
t100 = t96 * t93;
t99 = t96 * t95;
t92 = cos(pkin(6));
t84 = t92 * t100 + t101;
t90 = qJ(3) + pkin(11);
t88 = sin(t90);
t89 = cos(t90);
t98 = t88 * t103 - t84 * t89;
t97 = t89 * t103 + t84 * t88;
t86 = -t92 * t102 + t99;
t85 = t92 * t101 + t100;
t83 = t92 * t99 - t102;
t82 = t88 * t105 + t86 * t89;
t81 = t89 * t105 - t86 * t88;
t1 = [t98, -t85 * t89, t81, 0, 0, 0; t82, t83 * t89, -t97, 0, 0, 0; 0, t89 * t104, -t88 * t106 + t92 * t89, 0, 0, 0; t97, t85 * t88, -t82, 0, 0, 0; t81, -t83 * t88, t98, 0, 0, 0; 0, -t88 * t104, -t89 * t106 - t92 * t88, 0, 0, 0; t83, t86, 0, 0, 0, 0; t85, t84, 0, 0, 0, 0; 0, t106, 0, 0, 0, 0;];
JR_rot  = t1;

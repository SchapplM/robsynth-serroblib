% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRPRP5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2,theta4]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:46
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RPRPRP5_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP5_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRP5_jacobiR_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:45:53
% EndTime: 2019-02-26 20:45:53
% DurationCPUTime: 0.04s
% Computational Cost: add. (59->14), mult. (40->20), div. (0->0), fcn. (69->6), ass. (0->19)
t92 = pkin(9) + qJ(3);
t88 = sin(t92);
t93 = sin(qJ(1));
t100 = t93 * t88;
t91 = pkin(10) + qJ(5);
t89 = cos(t91);
t99 = t93 * t89;
t90 = cos(t92);
t98 = t93 * t90;
t94 = cos(qJ(1));
t97 = t94 * t88;
t96 = t94 * t89;
t95 = t94 * t90;
t87 = sin(t91);
t86 = t93 * t87 + t89 * t95;
t85 = t87 * t95 - t99;
t84 = -t94 * t87 + t89 * t98;
t83 = -t87 * t98 - t96;
t1 = [-t84, 0, -t88 * t96, 0, -t85, 0; t86, 0, -t88 * t99, 0, t83, 0; 0, 0, t90 * t89, 0, -t88 * t87, 0; -t100, 0, t95, 0, 0, 0; t97, 0, t98, 0, 0, 0; 0, 0, t88, 0, 0, 0; t83, 0, -t87 * t97, 0, t86, 0; t85, 0, -t87 * t100, 0, t84, 0; 0, 0, t90 * t87, 0, t88 * t89, 0;];
JR_rot  = t1;

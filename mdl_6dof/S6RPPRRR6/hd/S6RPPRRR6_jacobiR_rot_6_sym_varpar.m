% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPPRRR6
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:37
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RPPRRR6_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR6_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRR6_jacobiR_rot_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:37:35
% EndTime: 2019-02-26 20:37:36
% DurationCPUTime: 0.04s
% Computational Cost: add. (54->17), mult. (54->20), div. (0->0), fcn. (93->6), ass. (0->17)
t85 = sin(qJ(4));
t86 = sin(qJ(1));
t93 = t86 * t85;
t84 = qJ(5) + qJ(6);
t82 = sin(t84);
t87 = cos(qJ(4));
t92 = t87 * t82;
t83 = cos(t84);
t91 = t87 * t83;
t88 = cos(qJ(1));
t90 = t88 * t85;
t89 = t88 * t87;
t79 = -t86 * t82 + t83 * t90;
t78 = -t82 * t90 - t86 * t83;
t77 = -t88 * t82 - t83 * t93;
t76 = t82 * t93 - t88 * t83;
t1 = [t77, 0, 0, t83 * t89, t78, t78; t79, 0, 0, t86 * t91, -t76, -t76; 0, 0, 0, -t85 * t83, -t92, -t92; t76, 0, 0, -t82 * t89, -t79, -t79; t78, 0, 0, -t86 * t92, t77, t77; 0, 0, 0, t85 * t82, -t91, -t91; t86 * t87, 0, 0, t90, 0, 0; -t89, 0, 0, t93, 0, 0; 0, 0, 0, t87, 0, 0;];
JR_rot  = t1;

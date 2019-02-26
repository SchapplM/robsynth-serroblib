% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPPPR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta3,theta5]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:22
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRPPPR2_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR2_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPPR2_jacobiR_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:22:26
% EndTime: 2019-02-26 21:22:26
% DurationCPUTime: 0.03s
% Computational Cost: add. (57->12), mult. (40->20), div. (0->0), fcn. (69->6), ass. (0->19)
t83 = qJ(2) + pkin(9);
t79 = sin(t83);
t84 = sin(qJ(1));
t91 = t84 * t79;
t82 = pkin(10) + qJ(6);
t80 = cos(t82);
t90 = t84 * t80;
t81 = cos(t83);
t89 = t84 * t81;
t85 = cos(qJ(1));
t88 = t85 * t79;
t87 = t85 * t80;
t86 = t85 * t81;
t78 = sin(t82);
t77 = -t78 * t91 + t87;
t76 = t85 * t78 + t79 * t90;
t75 = t78 * t88 + t90;
t74 = -t84 * t78 + t79 * t87;
t1 = [t77, t78 * t86, 0, 0, 0, t74; t75, t78 * t89, 0, 0, 0, t76; 0, t79 * t78, 0, 0, 0, -t81 * t80; -t76, t80 * t86, 0, 0, 0, -t75; t74, t80 * t89, 0, 0, 0, t77; 0, t79 * t80, 0, 0, 0, t81 * t78; -t89, -t88, 0, 0, 0, 0; t86, -t91, 0, 0, 0, 0; 0, t81, 0, 0, 0, 0;];
JR_rot  = t1;

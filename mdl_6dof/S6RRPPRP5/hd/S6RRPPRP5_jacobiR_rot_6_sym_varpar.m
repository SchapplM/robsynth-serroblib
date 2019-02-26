% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPPRP5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta4]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:27
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRPPRP5_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP5_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPRP5_jacobiR_rot_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:27:27
% EndTime: 2019-02-26 21:27:27
% DurationCPUTime: 0.04s
% Computational Cost: add. (40->15), mult. (40->20), div. (0->0), fcn. (69->6), ass. (0->17)
t82 = sin(qJ(2));
t83 = sin(qJ(1));
t90 = t83 * t82;
t81 = pkin(9) + qJ(5);
t79 = sin(t81);
t84 = cos(qJ(2));
t89 = t84 * t79;
t80 = cos(t81);
t88 = t84 * t80;
t85 = cos(qJ(1));
t87 = t85 * t82;
t86 = t85 * t84;
t78 = -t79 * t90 + t85 * t80;
t77 = t85 * t79 + t80 * t90;
t76 = t79 * t87 + t83 * t80;
t75 = t83 * t79 - t80 * t87;
t1 = [t78, t79 * t86, 0, 0, -t75, 0; t76, t83 * t89, 0, 0, t77, 0; 0, t82 * t79, 0, 0, -t88, 0; -t83 * t84, -t87, 0, 0, 0, 0; t86, -t90, 0, 0, 0, 0; 0, t84, 0, 0, 0, 0; t77, -t80 * t86, 0, 0, t76, 0; t75, -t83 * t88, 0, 0, -t78, 0; 0, -t82 * t80, 0, 0, -t89, 0;];
JR_rot  = t1;

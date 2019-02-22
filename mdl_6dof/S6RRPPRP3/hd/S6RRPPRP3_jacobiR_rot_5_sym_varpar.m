% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPPRP3
% Use Code from Maple symbolic Code Generation
%
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 11:10
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRPPRP3_jacobiR_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP3_jacobiR_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RRPPRP3_jacobiR_rot_5_sym_varpar: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 11:10:44
% EndTime: 2019-02-22 11:10:44
% DurationCPUTime: 0.07s
% Computational Cost: add. (14->12), mult. (40->20), div. (0->0), fcn. (69->6), ass. (0->18)
t66 = sin(qJ(2));
t67 = sin(qJ(1));
t77 = t67 * t66;
t68 = cos(qJ(5));
t76 = t67 * t68;
t65 = sin(qJ(5));
t69 = cos(qJ(2));
t75 = t69 * t65;
t74 = t69 * t68;
t70 = cos(qJ(1));
t73 = t70 * t66;
t72 = t70 * t68;
t71 = t70 * t69;
t64 = -t67 * t65 + t66 * t72;
t63 = -t65 * t73 - t76;
t62 = -t70 * t65 - t66 * t76;
t61 = t65 * t77 - t72;
t1 = [t62, t68 * t71, 0, 0, t63, 0; t64, t67 * t74, 0, 0, -t61, 0; 0, t66 * t68, 0, 0, t75, 0; t61, -t65 * t71, 0, 0, -t64, 0; t63, -t67 * t75, 0, 0, t62, 0; 0, -t66 * t65, 0, 0, t74, 0; -t67 * t69, -t73, 0, 0, 0, 0; t71, -t77, 0, 0, 0, 0; 0, t69, 0, 0, 0, 0;];
JR_rot  = t1;

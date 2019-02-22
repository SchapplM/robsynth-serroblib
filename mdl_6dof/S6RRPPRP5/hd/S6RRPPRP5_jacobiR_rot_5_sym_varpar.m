% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
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
% Datum: 2019-02-22 11:12
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRPPRP5_jacobiR_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP5_jacobiR_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPRP5_jacobiR_rot_5_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 11:12:01
% EndTime: 2019-02-22 11:12:01
% DurationCPUTime: 0.06s
% Computational Cost: add. (36->11), mult. (40->20), div. (0->0), fcn. (69->6), ass. (0->17)
t71 = sin(qJ(2));
t72 = sin(qJ(1));
t79 = t72 * t71;
t70 = pkin(9) + qJ(5);
t68 = sin(t70);
t73 = cos(qJ(2));
t78 = t73 * t68;
t69 = cos(t70);
t77 = t73 * t69;
t74 = cos(qJ(1));
t76 = t74 * t71;
t75 = t74 * t73;
t67 = -t68 * t79 + t74 * t69;
t66 = t74 * t68 + t69 * t79;
t65 = t68 * t76 + t72 * t69;
t64 = -t72 * t68 + t69 * t76;
t1 = [t67, t68 * t75, 0, 0, t64, 0; t65, t72 * t78, 0, 0, t66, 0; 0, t71 * t68, 0, 0, -t77, 0; -t66, t69 * t75, 0, 0, -t65, 0; t64, t72 * t77, 0, 0, t67, 0; 0, t71 * t69, 0, 0, t78, 0; -t72 * t73, -t76, 0, 0, 0, 0; t75, -t79, 0, 0, 0, 0; 0, t73, 0, 0, 0, 0;];
JR_rot  = t1;

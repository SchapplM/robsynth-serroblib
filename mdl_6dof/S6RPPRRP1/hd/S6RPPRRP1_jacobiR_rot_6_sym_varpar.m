% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPPRRP1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta2,theta3]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 10:13
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RPPRRP1_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP1_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRP1_jacobiR_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 10:13:25
% EndTime: 2019-02-22 10:13:25
% DurationCPUTime: 0.04s
% Computational Cost: add. (59->14), mult. (40->20), div. (0->0), fcn. (69->6), ass. (0->19)
t80 = qJ(1) + pkin(9);
t76 = sin(t80);
t81 = sin(qJ(5));
t88 = t76 * t81;
t82 = cos(qJ(5));
t87 = t76 * t82;
t79 = pkin(10) + qJ(4);
t77 = cos(t79);
t86 = t77 * t81;
t85 = t77 * t82;
t78 = cos(t80);
t84 = t78 * t81;
t83 = t78 * t82;
t75 = sin(t79);
t74 = t77 * t83 + t88;
t73 = -t77 * t84 + t87;
t72 = -t76 * t85 + t84;
t71 = t76 * t86 + t83;
t1 = [t72, 0, 0, -t75 * t83, t73, 0; t74, 0, 0, -t75 * t87, -t71, 0; 0, 0, 0, t85, -t75 * t81, 0; t71, 0, 0, t75 * t84, -t74, 0; t73, 0, 0, t75 * t88, t72, 0; 0, 0, 0, -t86, -t75 * t82, 0; -t76 * t75, 0, 0, t78 * t77, 0, 0; t78 * t75, 0, 0, t76 * t77, 0, 0; 0, 0, 0, t75, 0, 0;];
JR_rot  = t1;

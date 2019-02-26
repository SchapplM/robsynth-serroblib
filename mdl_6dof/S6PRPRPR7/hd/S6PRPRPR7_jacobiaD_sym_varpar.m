% Zeitableitung der analytischen Jacobi-Matrix für beliebiges Segment von
% S6PRPRPR7
% Use Code from Maple symbolic Code Generation
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% Zeitableitung der Winkeldarstellung des Endeffektors in Basis-Koordinaten
% 
% Winkeldarstellung: Euler-XYZ-Winkel, rotx(alpha)*roty(beta)*rotz(gamma)
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt.
%   Wie in S6PRPRPR7_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1]';
% 
% Output:
% JaD [6x6]
%   Zeitableitung der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 14:59
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function JaD = S6PRPRPR7_jacobiaD_floatb_twist_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)



%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR7_jacobiaD_floatb_twist_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR7_jacobiaD_floatb_twist_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRPRPR7_jacobiaD_floatb_twist_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRPRPR7_jacobiaD_floatb_twist_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRPRPR7_jacobiaD_floatb_twist_sym_varpar: pkin has to be [10x1] (double)');
%% Function calls
if link_index == 0
	JaD=S6PRPRPR7_jacobiaD_0_floatb_twist_sym_varpar(qJ, qJD, r_i_i_C, pkin);
elseif link_index == 1
	JaD=S6PRPRPR7_jacobiaD_1_floatb_twist_sym_varpar(qJ, qJD, r_i_i_C, pkin);
elseif link_index == 2
	JaD=S6PRPRPR7_jacobiaD_2_floatb_twist_sym_varpar(qJ, qJD, r_i_i_C, pkin);
elseif link_index == 3
	JaD=S6PRPRPR7_jacobiaD_3_floatb_twist_sym_varpar(qJ, qJD, r_i_i_C, pkin);
elseif link_index == 4
	JaD=S6PRPRPR7_jacobiaD_4_floatb_twist_sym_varpar(qJ, qJD, r_i_i_C, pkin);
elseif link_index == 5
	JaD=S6PRPRPR7_jacobiaD_5_floatb_twist_sym_varpar(qJ, qJD, r_i_i_C, pkin);
elseif link_index == 6
	JaD=S6PRPRPR7_jacobiaD_6_floatb_twist_sym_varpar(qJ, qJD, r_i_i_C, pkin);
else
	JaD=NaN(6,6);
end